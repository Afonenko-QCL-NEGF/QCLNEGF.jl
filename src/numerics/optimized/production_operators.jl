# Production finite open-chain recursion. Each logical worker owns its scratch;
# task migration cannot alias another worker's matrices. The reference oracle
# in reference/periodicity.jl remains independent and unchanged.
"""Owned LAPACK LU/inverse buffers, reused for independent small dense blocks."""
struct _ProductionInverseWorkspace
    pivots::Vector{LinearAlgebra.BlasInt}
    work::Vector{ComplexF64}
    info::Vector{LinearAlgebra.BlasInt}
end
_ProductionInverseWorkspace(n::Int) = _ProductionInverseWorkspace(
    zeros(LinearAlgebra.BlasInt, n),
    zeros(ComplexF64, max(1, 64n)),
    zeros(LinearAlgebra.BlasInt, 1),
)

function _production_inverse!(A::Matrix{ComplexF64}, workspace::_ProductionInverseWorkspace)
    n = LinearAlgebra.BlasInt(size(A, 1))
    size(A, 2) == n == length(workspace.pivots) ||
        throw(DimensionMismatch("inverse workspace differs from matrix"))
    ccall(
        (LinearAlgebra.BLAS.@blasfunc(zgetrf_), LinearAlgebra.LAPACK.libblastrampoline),
        Cvoid,
        (
            Ref{LinearAlgebra.BlasInt},
            Ref{LinearAlgebra.BlasInt},
            Ptr{ComplexF64},
            Ref{LinearAlgebra.BlasInt},
            Ptr{LinearAlgebra.BlasInt},
            Ptr{LinearAlgebra.BlasInt},
        ),
        n,
        n,
        A,
        n,
        workspace.pivots,
        workspace.info,
    )
    info = workspace.info[1]
    info > 0 && throw(SingularException(info))
    info < 0 && throw(ArgumentError("LAPACK getrf rejected argument $(-info)"))
    lwork = LinearAlgebra.BlasInt(length(workspace.work))
    ccall(
        (LinearAlgebra.BLAS.@blasfunc(zgetri_), LinearAlgebra.LAPACK.libblastrampoline),
        Cvoid,
        (
            Ref{LinearAlgebra.BlasInt},
            Ptr{ComplexF64},
            Ref{LinearAlgebra.BlasInt},
            Ptr{LinearAlgebra.BlasInt},
            Ptr{ComplexF64},
            Ref{LinearAlgebra.BlasInt},
            Ptr{LinearAlgebra.BlasInt},
        ),
        n,
        A,
        n,
        workspace.pivots,
        workspace.work,
        lwork,
        workspace.info,
    )
    info = workspace.info[1]
    info > 0 && throw(SingularException(info))
    info < 0 && throw(ArgumentError("LAPACK getri rejected argument $(-info)"))
    return A
end

# A partial-pivot LU and triangular solves for tiny blocks avoid millions of
# dispatches into general LAPACK kernels. This is not a closed-form inverse:
# pivots, singular checks and ordinary Float64 arithmetic remain explicit.
# Only the cavity recursion selects this path. Seed/Dyson retain their shared
# LAPACK arithmetic and bitwise reference contract. Larger bases also retain
# LAPACK. The caller owns all scratch storage.
function _cavity_inverse!(A::Matrix{ComplexF64}, workspace::_ProductionInverseWorkspace)
    n = size(A, 1)
    size(A, 2) == n == length(workspace.pivots) ||
        throw(DimensionMismatch("inverse workspace differs from matrix"))
    1 <= n <= 8 || return _production_inverse!(A, workspace)
    # The common five-state basis has a statically sized scalar kernel.
    n == 5 && return _production_inverse_small!(A, workspace, Val(5))
    return _production_inverse_small!(A, workspace, Val(n))
end

function _production_inverse_small!(A, workspace, ::Val{n}) where {n}
    @inbounds for k = 1:n
        pivot = k
        largest = abs(real(A[k, k])) + abs(imag(A[k, k]))
        for i = (k+1):n
            magnitude = abs(real(A[i, k])) + abs(imag(A[i, k]))
            if magnitude > largest
                largest, pivot = magnitude, i
            end
        end
        iszero(A[pivot, k]) && throw(SingularException(k))
        workspace.pivots[k] = pivot
        if pivot != k
            for j = 1:n
                A[k, j], A[pivot, j] = A[pivot, j], A[k, j]
            end
        end
        for i = (k+1):n
            A[i, k] /= A[k, k]
        end
        for j = (k+1):n, i = (k+1):n
            A[i, j] -= A[i, k] * A[k, j]
        end
    end
    # Solve A_original * X = I in the preallocated LAPACK work vector.
    # Apply the same row swaps to each identity column, then solve L and U.
    X = workspace.work
    @inbounds for column = 1:n
        offset = (column-1)*n
        for i = 1:n
            X[offset+i] = i == column ? 1.0 : 0.0
        end
        for k = 1:n
            pivot = workspace.pivots[k]
            X[offset+k], X[offset+pivot] = X[offset+pivot], X[offset+k]
        end
        for i = 1:n
            value = X[offset+i]
            for k = 1:(i-1)
                value -= A[i, k] * X[offset+k]
            end
            X[offset+i] = value
        end
        for i = n:-1:1
            value = X[offset+i]
            for k = (i+1):n
                value -= A[i, k] * X[offset+k]
            end
            X[offset+i] = value / A[i, i]
        end
    end
    @inbounds for j = 1:n, i = 1:n
        A[i, j] = X[i+(j-1)*n]
    end
    return A
end

struct _CavityWorkspace
    retarded::Matrix{ComplexF64}
    lesser::Matrix{ComplexF64}
    greater::Matrix{ComplexF64}
    denominator::Matrix{ComplexF64}
    local_lesser::Matrix{ComplexF64}
    local_greater::Matrix{ComplexF64}
    temporary::Matrix{ComplexF64}
    product::Matrix{ComplexF64}
    inverse::_ProductionInverseWorkspace
end
_CavityWorkspace(n::Int) = _CavityWorkspace(
    (zeros(ComplexF64, n, n) for _ = 1:8)...,
    _ProductionInverseWorkspace(n),
)

# Compact precomputed remaps also cover nonuniform conservative quadrature.
function _cavity_remaps(problem, periods, energy_shift)
    nodes, weights = problem.grids.ε, problem.grids.wᴱ
    count = length(nodes)
    drop = _scaled_physics(problem.physical, problem.scales).Eᵖ
    edges = vcat(first(nodes), (nodes[1:(end-1)] .+ nodes[2:end]) ./ 2, last(nodes))
    offsets = collect((-periods):periods)
    return map(offsets) do displacement_index
        displacement = displacement_index * drop
        rows = Vector{Vector{Pair{Int,Float64}}}(undef, count)
        for e = 1:count
            row = Pair{Int,Float64}[]
            if energy_shift === :conservative_pair
                left, right = edges[e] + displacement, edges[e+1] + displacement
                lo = clamp(searchsortedlast(edges, left), 1, count)
                hi = clamp(searchsortedfirst(edges, right), 1, count)
                for source = lo:hi
                    overlap =
                        max(0.0, min(right, edges[source+1]) - max(left, edges[source]))
                    overlap > 0 && push!(row, source => overlap / weights[e])
                end
            else
                target = nodes[e] + displacement
                if first(nodes) <= target <= last(nodes)
                    if target == last(nodes)
                        push!(row, count => 1.0)
                    else
                        j = clamp(searchsortedlast(nodes, target), 1, count-1)
                        fraction = (target - nodes[j]) / (nodes[j+1] - nodes[j])
                        push!(row, j => 1-fraction, j+1 => fraction)
                    end
                end
            end
            rows[e] = row
        end
        rows
    end
end

function _cavity_local!(output, field, row, m)
    fill!(output, 0)
    @inbounds for (source, weight) in row, b in axes(output, 2), a in axes(output, 1)
        output[a, b] += weight * field[source, m, a, b]
    end
    return output
end

function _cavity_sandwich!(output, coupling, matrix, temporary)
    if size(output, 1) == 5
        return _cavity_sandwich_small!(output, coupling, matrix, temporary, Val(5))
    elseif size(output, 1) <= 8
        _production_matrix_product!(temporary, coupling, matrix)
        _production_matrix_product_adjoint_right!(output, temporary, coupling)
    else
        mul!(temporary, coupling, matrix)
        mul!(output, temporary, adjoint(coupling))
    end
    return output
end

function _cavity_sandwich_small!(output, coupling, matrix, temporary, ::Val{n}) where {n}
    @inbounds for j = 1:n, i = 1:n
        value = 0.0 + 0.0im
        for k = 1:n
            value += coupling[i, k] * matrix[k, j]
        end
        temporary[i, j] = value
    end
    @inbounds for j = 1:n, i = 1:n
        value = 0.0 + 0.0im
        for k = 1:n
            value += temporary[i, k] * conj(coupling[j, k])
        end
        output[i, j] = value
    end
    return output
end

function _cavity_side!(
    w::_CavityWorkspace,
    family,
    h,
    remaps,
    e,
    m,
    nodes,
    drop,
    periods,
    coupling,
    direction,
    destination,
)
    # The outside of the finite chain has exactly zero embedding; skip its
    # three zero matrix sandwiches. Later cells use the computed cavity state.
    boundary = true
    indices = direction === :minus ? ((-periods):-1) : (periods:-1:1)
    for j in indices
        row = remaps[j+periods+1][e]
        _cavity_local!(w.denominator, family.Σᴿ, row, m)
        boundary || _cavity_sandwich!(w.product, coupling, w.retarded, w.temporary)
        @inbounds for b in axes(w.denominator, 2), a in axes(w.denominator, 1)
            w.denominator[a, b] =
                -h[m, a, b] - w.denominator[a, b] -
                (boundary ? zero(ComplexF64) : w.product[a, b])
        end
        for a in axes(w.denominator, 1)
            w.denominator[a, a] += nodes[e] + j*drop
        end
        _cavity_local!(w.local_lesser, family.Σˡ, row, m)
        if !boundary
            _cavity_sandwich!(w.product, coupling, w.lesser, w.temporary)
            w.local_lesser .+= w.product
        end
        _cavity_local!(w.local_greater, family.Σᵍ, row, m)
        if !boundary
            _cavity_sandwich!(w.product, coupling, w.greater, w.temporary)
            w.local_greater .+= w.product
        end
        # LU works on bounded Nb×Nb storage. No central-cell inverse is needed
        # when the caller requests only the two embeddings.
        copyto!(w.retarded, _cavity_inverse!(w.denominator, w.inverse))
        _cavity_sandwich!(w.lesser, w.retarded, w.local_lesser, w.temporary)
        _cavity_sandwich!(w.greater, w.retarded, w.local_greater, w.temporary)
        boundary = false
    end
    for (field, matrix) in (
        (destination.Σᴿ, w.retarded),
        (destination.Σˡ, w.lesser),
        (destination.Σᵍ, w.greater),
    )
        _cavity_sandwich!(w.product, coupling, matrix, w.temporary)
        for b in axes(w.product, 2), a in axes(w.product, 1)
            field[e, m, a, b] = w.product[a, b]
        end
    end
    return nothing
end

"""Immutable geometry/Hartree data scoped to a single SCBA solve.

No self-energy or iteration output is retained. The original Hamiltonian is
copied so accidental in-place changes fail closed instead of reusing stale data.
"""
function _cavity_production_plan(problem, h, options)
    periods = options.algorithms.embedding_periods
    periods >= 1 || throw(ArgumentError("cavity requires neighbouring cells"))
    energy_shift = options.algorithms.energy_shift
    energy_shift in (:dense, :sparse_plan, :conservative_pair) ||
        throw(ArgumentError("unsupported cavity energy shift"))
    NE, Nk, Nb = problem.numerical.N_E, problem.numerical.N_k, problem.numerical.N_b
    size(h) == (Nk, Nb, Nb) || throw(DimensionMismatch("cavity Hamiltonian differs"))
    length(problem.grids.ε) == NE || throw(DimensionMismatch("cavity energy grid differs"))
    hlocal = Array{ComplexF64}(undef, size(h))
    for m = 1:Nk
        block = Matrix(view(h, m, :, :))
        all(isfinite, block) &&
        norm(block-block') <= 128eps(Float64)*max(norm(block), 1.0) ||
            throw(ArgumentError("onsite blocks must be Hermitian"))
        hlocal[m, :, :] .= (block+block')/2
    end
    return (;
        source_grids = problem.grids,
        source_energy = copy(problem.grids.ε),
        source_weights = copy(problem.grids.wᴱ),
        source_h = copy(h),
        hlocal,
        remaps = _cavity_remaps(problem, periods, energy_shift),
        drop = _scaled_physics(problem.physical, problem.scales).Eᵖ,
        periods,
        energy_shift,
    )
end

function _cavity_embedding_production(
    problem::NEGFProblem,
    family::SelfEnergyFamily,
    h::AbstractArray{<:Number,3},
    options::ProductionOptions;
    plan = nothing,
)
    plan === nothing && (plan = _cavity_production_plan(problem, h, options))
    NE, Nk, Nb, Nb2 = size(family.Σᴿ)
    (NE, Nk, Nb) == (problem.numerical.N_E, problem.numerical.N_k, problem.numerical.N_b) &&
    Nb == Nb2 || throw(DimensionMismatch("cavity dimensions differ from problem"))
    size(family.Σˡ) == size(family.Σᵍ) == size(family.Σᴿ) ||
        throw(DimensionMismatch("cavity components differ"))
    plan.source_grids === problem.grids &&
    plan.source_h == h &&
    plan.source_energy == problem.grids.ε &&
    plan.source_weights == problem.grids.wᴱ &&
    plan.periods == options.algorithms.embedding_periods &&
    plan.energy_shift === options.algorithms.energy_shift &&
    plan.drop == _scaled_physics(problem.physical, problem.scales).Eᵖ ||
        throw(ArgumentError("cavity plan is stale; rebuild it for this SCBA solve"))
    (; hlocal, remaps, periods, drop) = plan
    plus, minus = _empty_selfenergy(NE, Nk, Nb), _empty_selfenergy(NE, Nk, Nb)
    workers = _production_worker_count(options, NE*Nk)
    function evaluate_worker(worker)
        workspace = _CavityWorkspace(Nb)
        first_linear = fld((worker - 1) * NE * Nk, workers) + 1
        last_linear = fld(worker * NE * Nk, workers)
        for linear = first_linear:last_linear
            e, m = (linear-1) % NE + 1, (linear-1) ÷ NE + 1
            _cavity_side!(
                workspace,
                family,
                hlocal,
                remaps,
                e,
                m,
                plan.source_energy,
                drop,
                periods,
                adjoint(problem.basis.T₊),
                :minus,
                minus,
            )
            _cavity_side!(
                workspace,
                family,
                hlocal,
                remaps,
                e,
                m,
                plan.source_energy,
                drop,
                periods,
                problem.basis.T₊,
                :plus,
                plus,
            )
        end
    end
    if workers == 1
        evaluate_worker(1)
    else
        @sync for worker = 1:workers
            Base.Threads.@spawn evaluate_worker(worker)
        end
    end
    total = SelfEnergyFamily(plus.Σᴿ+minus.Σᴿ, plus.Σˡ+minus.Σˡ, plus.Σᵍ+minus.Σᵍ)
    return total, plus, minus
end

function _production_static_contraction!(
    Σ::AbstractArray{ComplexF64,4},
    op::DenseProductionKernel,
    G::AbstractArray{<:Number,4},
    wᵏ::AbstractVector{<:Real},
    qᴷ::Real,
    energy_chunk::Integer;
    α::Number = 1,
    β::Number = 0,
    parallel_backend::Symbol = :blas,
    worker_count::Integer = 0,
)
    NE, Nk, Nb, Nb2 = size(G)
    (Nk, Nb) == (op.N_k, op.N_b) && Nb == Nb2 ||
        throw(DimensionMismatch("Green field and production kernel differ"))
    size(Σ) == size(G) || throw(DimensionMismatch("output shape differs"))
    length(wᵏ) == Nk || throw(DimensionMismatch("radial weights differ"))
    energy_chunk > 0 || throw(ArgumentError("energy_chunk must be positive"))
    Ncompound = Nk * Nb^2
    G2 = reshape(G, NE, Ncompound)
    Σ2 = reshape(Σ, NE, Ncompound)
    weights = repeat(Float64.(wᵏ), Nb^2)
    first_energies = collect(1:energy_chunk:NE)
    workers =
        _production_worker_count(parallel_backend, worker_count, length(first_energies))
    chunk_capacity = min(Int(energy_chunk), NE)
    factor = α * qᴷ

    function run_worker(worker::Int)
        X = Matrix{ComplexF64}(undef, Ncompound, chunk_capacity)
        Y = similar(X)
        for chunk_index = worker:workers:length(first_energies)
            first_energy = first_energies[chunk_index]
            last_energy = min(first_energy + energy_chunk - 1, NE)
            energies = first_energy:last_energy
            Nc = length(energies)
            Xview = @view X[:, 1:Nc]
            Yview = @view Y[:, 1:Nc]
            @inbounds for j = 1:Nc, q = 1:Ncompound
                Xview[q, j] = weights[q] * G2[energies[j], q]
            end
            mul!(Yview, op.matrix, Xview)
            if iszero(β)
                @inbounds for j = 1:Nc, q = 1:Ncompound
                    Σ2[energies[j], q] = factor * Yview[q, j]
                end
            else
                @inbounds for j = 1:Nc, q = 1:Ncompound
                    Σ2[energies[j], q] = β * Σ2[energies[j], q] + factor * Yview[q, j]
                end
            end
        end
        return nothing
    end

    if workers == 1
        run_worker(1)
    else
        Base.Threads.@threads :static for worker = 1:workers
            run_worker(worker)
        end
    end
    return Σ
end

function _production_static_contraction!(
    Σ::AbstractArray{ComplexF64,4},
    op::MomentumIndependentProductionKernel,
    G::AbstractArray{<:Number,4},
    wᵏ::AbstractVector{<:Real},
    qᴷ::Real,
    energy_chunk::Integer;
    α::Number = 1,
    β::Number = 0,
    parallel_backend::Symbol = :blas,
    worker_count::Integer = 0,
)
    NE, Nk, Nb, Nb2 = size(G)
    (Nk, Nb) == (op.N_k, op.N_b) && Nb == Nb2 ||
        throw(DimensionMismatch("Green field and production kernel differ"))
    size(Σ) == size(G) || throw(DimensionMismatch("output shape differs"))
    length(wᵏ) == Nk || throw(DimensionMismatch("radial weights differ"))
    energy_chunk > 0 || throw(ArgumentError("energy_chunk must be positive"))
    Nstate = Nb^2
    first_energies = collect(1:energy_chunk:NE)
    workers =
        _production_worker_count(parallel_backend, worker_count, length(first_energies))
    chunk_capacity = min(Int(energy_chunk), NE)
    factor = α * qᴷ

    function run_worker(worker::Int)
        X = zeros(ComplexF64, Nstate, chunk_capacity)
        Y = similar(X)
        for chunk_index = worker:workers:length(first_energies)
            first_energy = first_energies[chunk_index]
            last_energy = min(first_energy + energy_chunk - 1, NE)
            energies = first_energy:last_energy
            Nc = length(energies)
            Xview = @view X[:, 1:Nc]
            Yview = @view Y[:, 1:Nc]
            fill!(Xview, 0)
            @inbounds for j = 1:Nc, d = 1:Nb, c = 1:Nb, mp = 1:Nk
                q = c + Nb * (d - 1)
                Xview[q, j] += wᵏ[mp] * G[energies[j], mp, c, d]
            end
            mul!(Yview, op.block, Xview)
            @inbounds for j = 1:Nc, b = 1:Nb, a = 1:Nb, m = 1:Nk
                q = a + Nb * (b - 1)
                e = energies[j]
                Σ[e, m, a, b] =
                    iszero(β) ? factor * Yview[q, j] :
                    β * Σ[e, m, a, b] + factor * Yview[q, j]
            end
        end
        return nothing
    end

    if workers == 1
        run_worker(1)
    else
        Base.Threads.@threads :static for worker = 1:workers
            run_worker(worker)
        end
    end
    return Σ
end

function _production_static_contraction!(
    Σ::AbstractArray{ComplexF64,4},
    op::LowRankProductionKernel,
    G::AbstractArray{<:Number,4},
    wᵏ::AbstractVector{<:Real},
    qᴷ::Real,
    energy_chunk::Integer;
    α::Number = 1,
    β::Number = 0,
    parallel_backend::Symbol = :blas,
    worker_count::Integer = 0,
)
    NE, Nk, Nb, Nb2 = size(G)
    (Nk, Nb) == (op.N_k, op.N_b) && Nb == Nb2 ||
        throw(DimensionMismatch("Green field and low-rank kernel differ"))
    size(Σ) == size(G) || throw(DimensionMismatch("output shape differs"))
    length(wᵏ) == Nk || throw(DimensionMismatch("radial weights differ"))
    energy_chunk > 0 || throw(ArgumentError("energy_chunk must be positive"))
    Ncompound = Nk * Nb^2
    G2 = reshape(G, NE, Ncompound)
    Σ2 = reshape(Σ, NE, Ncompound)
    weights = repeat(Float64.(wᵏ), Nb^2)
    first_energies = collect(1:energy_chunk:NE)
    workers =
        _production_worker_count(parallel_backend, worker_count, length(first_energies))
    capacity = min(Int(energy_chunk), NE)
    factor = α * qᴷ

    function run_worker(worker::Int)
        X = Matrix{ComplexF64}(undef, Ncompound, capacity)
        reduced = Matrix{ComplexF64}(undef, op.rank, capacity)
        Y = Matrix{ComplexF64}(undef, Ncompound, capacity)
        for chunk_index = worker:workers:length(first_energies)
            first_energy = first_energies[chunk_index]
            last_energy = min(first_energy + energy_chunk - 1, NE)
            energies = first_energy:last_energy
            count = length(energies)
            Xview = @view X[:, 1:count]
            reduced_view = @view reduced[:, 1:count]
            Yview = @view Y[:, 1:count]
            @inbounds for j = 1:count, q = 1:Ncompound
                Xview[q, j] = weights[q] * G2[energies[j], q]
            end
            mul!(reduced_view, op.right_adjoint, Xview)
            mul!(Yview, op.left, reduced_view)
            @inbounds for j = 1:count, q = 1:Ncompound
                e = energies[j]
                Σ2[e, q] =
                    iszero(β) ? factor * Yview[q, j] : β * Σ2[e, q] + factor * Yview[q, j]
            end
        end
        return nothing
    end
    if workers == 1
        run_worker(1)
    else
        Base.Threads.@threads :static for worker = 1:workers
            run_worker(worker)
        end
    end
    return Σ
end

function _production_static_contraction!(
    ::AbstractArray{ComplexF64,4},
    ::LiteralProductionKernel,
    ::AbstractArray{<:Number,4},
    ::AbstractVector{<:Real},
    ::Real,
    ::Integer;
    kwargs...,
)
    throw(ArgumentError("a LiteralProductionKernel requires its source six-index tensor"))
end

"""
    production_static_contraction(operator, G, wᵏ, qᴷ=1; energy_chunk=128)

Memory-bounded BLAS evaluation of the *same* six-index contraction as
`_static_contraction`.  The dense operator changes storage order only.  A
momentum-independent kernel is reduced exactly before multiplication.

See [Microscopic kernels](@ref theory-kernels) and the
[production BLAS contraction](@ref theory-production).
"""
function production_static_contraction(
    op::AbstractProductionKernel,
    G::AbstractArray{<:Number,4},
    wᵏ::AbstractVector{<:Real},
    qᴷ::Real = 1;
    energy_chunk::Integer = 128,
    parallel_backend::Symbol = :blas,
    worker_count::Integer = 0,
)
    Σ = zeros(ComplexF64, size(G))
    return _production_static_contraction!(
        Σ,
        op,
        G,
        wᵏ,
        qᴷ,
        energy_chunk;
        parallel_backend,
        worker_count,
    )
end

function production_static_contraction(
    K::AbstractArray{<:Number,6},
    G::AbstractArray{<:Number,4},
    wᵏ::AbstractVector{<:Real},
    qᴷ::Real = 1;
    energy_chunk::Integer = 128,
    parallel_backend::Symbol = :blas,
    worker_count::Integer = 0,
)
    return production_static_contraction(
        _production_kernel(K),
        G,
        wᵏ,
        qᴷ;
        energy_chunk,
        parallel_backend,
        worker_count,
    )
end

function _uniform_energy_spacing(ε::AbstractVector{<:Real})
    length(ε) ≥ 2 || throw(ArgumentError("at least two energy nodes are required"))
    Δε = Float64(ε[2] - ε[1])
    Δε > 0 || throw(ArgumentError("energy grid must be strictly increasing"))
    deviation = maximum(abs.(diff(ε) .- Δε))
    tolerance = 256eps(Float64) * max(maximum(abs, ε), abs(Δε), 1.0)
    deviation ≤ tolerance ||
        throw(ArgumentError("FFT Hilbert transform requires a uniform energy grid"))
    return Δε
end

function _fft_hilbert_with_length(
    Γ::AbstractArray{<:Number,4},
    ε::AbstractVector{<:Real},
    wᴱ::AbstractVector{<:Real},
    L::Integer,
    column_chunk::Integer,
)
    NE, Nk, Nb, Nb2 = size(Γ)
    Nb == Nb2 || throw(DimensionMismatch("Γ blocks must be square"))
    length(ε) == NE == length(wᴱ) || throw(DimensionMismatch("energy data differ"))
    L ≥ 2NE - 1 || throw(ArgumentError("FFT padding is too short"))
    Δε = _uniform_energy_spacing(ε)
    kernel = zeros(ComplexF64, L)
    for lag = (-(NE-1)):(NE-1)
        iszero(lag) && continue
        kernel[lag+NE] = inv(2π * lag * Δε)
    end
    FFTW.fft!(kernel)
    Ncolumns = Nk * Nb^2
    Γ2 = reshape(Γ, NE, Ncolumns)
    Λ2 = Matrix{ComplexF64}(undef, NE, Ncolumns)
    for first_column = 1:column_chunk:Ncolumns
        last_column = min(first_column + column_chunk - 1, Ncolumns)
        columns = first_column:last_column
        Nc = length(columns)
        buffer = zeros(ComplexF64, L, Nc)
        for j = 1:Nc, e = 1:NE
            buffer[e, j] = wᴱ[e] * Γ2[e, columns[j]]
        end
        FFTW.fft!(buffer, 1)
        for j = 1:Nc, i = 1:L
            buffer[i, j] *= kernel[i]
        end
        FFTW.ifft!(buffer, 1)
        for j = 1:Nc, e = 1:NE
            Λ2[e, columns[j]] = buffer[e+NE-1, j]
        end
    end
    return reshape(Λ2, size(Γ))
end

"""
    fft_hilbert_transform(Γ, ε, wᴱ; column_chunk=32,
                          return_roundoff=false)

`O(N_E log N_E)` zero-padded *linear* convolution for the exact discrete
principal-value sum used by [`direct_hilbert_transform`](@ref QCLNEGF.QCLReferenceOperators.direct_hilbert_transform).  It is not a
cyclic Hilbert transform and introduces no edge wrap.  The only extra
difference from the direct formula is floating-point summation order.

With `return_roundoff=true`, a second valid transform with twice the padding
length provides an empirical rounding-order residual.

See [the direct Hilbert relation](@ref eq-direct-hilbert),
[Microscopic kernels](@ref theory-kernels), and
[Production backend](@ref theory-production).
"""
function fft_hilbert_transform(
    Γ::AbstractArray{<:Number,4},
    ε::AbstractVector{<:Real},
    wᴱ::AbstractVector{<:Real};
    column_chunk::Integer = 32,
    return_roundoff::Bool = false,
)
    column_chunk > 0 || throw(ArgumentError("column_chunk must be positive"))
    NE = size(Γ, 1)
    L = nextpow(2, 2NE - 1)
    Λ = _fft_hilbert_with_length(Γ, ε, wᴱ, L, column_chunk)
    return_roundoff || return Λ
    Λsecond = _fft_hilbert_with_length(Γ, ε, wᴱ, 2L, column_chunk)
    Λsecond .-= Λ
    residual = norm(Λsecond) / (norm(Λ) + 1e-14)
    return Λ, residual
end

"""
Retarded self-energy reconstructed with the exact-discrete FFT Hilbert map.

See [the direct Hilbert relation](@ref eq-direct-hilbert),
[Microscopic kernels](@ref theory-kernels), and
[Production backend](@ref theory-production).
"""
function retarded_self_energy_fft(
    Σˡ::AbstractArray{<:Number,4},
    Σᵍ::AbstractArray{<:Number,4},
    grids::ModelGrids;
    column_chunk::Integer = 32,
    return_roundoff::Bool = false,
)
    size(Σˡ) == size(Σᵍ) || throw(DimensionMismatch("self-energy components differ"))
    Γ = im .* (Σᵍ .- Σˡ)
    if return_roundoff
        Λ, residual = fft_hilbert_transform(
            Γ,
            grids.ε,
            grids.wᴱ;
            column_chunk,
            return_roundoff = true,
        )
        return ComplexF64.(Λ .- 0.5im .* Γ), ComplexF64.(Γ), residual
    end
    Λ = fft_hilbert_transform(Γ, grids.ε, grids.wᴱ; column_chunk)
    return ComplexF64.(Λ .- 0.5im .* Γ), ComplexF64.(Γ)
end

function _retarded_self_energy_fft_production(
    Σˡ::AbstractArray{<:Number,4},
    Σᵍ::AbstractArray{<:Number,4},
    grids::ModelGrids,
    plan,
    options::ProductionOptions;
    return_roundoff::Bool = false,
)
    size(Σˡ) == size(Σᵍ) || throw(DimensionMismatch("self-energy components differ"))
    Γ = Array{ComplexF64}(undef, size(Σˡ))
    _production_parallel_linear!(length(Γ), options) do index
        Γ[index] = im * (Σᵍ[index] - Σˡ[index])
    end
    if return_roundoff
        Λ, residual = production_fft_hilbert_transform(
            Γ,
            grids.ε,
            grids.wᴱ,
            plan;
            return_roundoff = true,
            worker_count = options.worker_count,
        )
        # The transform output is exclusively owned by this call; after the
        # roundoff comparison its real-part role is finished. Reuse it as ΣR.
        Σᴿ = Λ
        _production_parallel_linear!(length(Γ), options) do index
            Σᴿ[index] = Λ[index] - 0.5im * Γ[index]
        end
        return Σᴿ, Γ, residual
    end
    Λ = production_fft_hilbert_transform(
        Γ,
        grids.ε,
        grids.wᴱ,
        plan;
        worker_count = options.worker_count,
    )
    Σᴿ = Λ
    _production_parallel_linear!(length(Γ), options) do index
        Σᴿ[index] = Λ[index] - 0.5im * Γ[index]
    end
    return Σᴿ, Γ
end

function _embedding_component_production(
    T::AbstractMatrix,
    G::AbstractArray{<:Number,4},
    plan::EnergyShiftPlan,
    options::ProductionOptions = ProductionOptions(),
)
    NE, Nk, Nb, Nb2 = size(G)
    Nb == Nb2 && size(T) == (Nb, Nb) ||
        throw(DimensionMismatch("embedding dimensions differ"))
    Σ = zeros(ComplexF64, size(G))
    blocks = NE * Nk
    workers = _production_worker_count(options, blocks)
    _production_worker_ranges!(blocks, workers) do _, first_linear, last_linear
        local shifted = zeros(ComplexF64, Nb, Nb)
        local temporary = similar(shifted)
        local result = similar(shifted)
        @inbounds for linear = first_linear:last_linear
            local e = (linear - 1) % NE + 1
            local m = (linear - 1) ÷ NE + 1
            plan.valid[e] || continue
            local l = plan.left[e]
            local r = plan.right[e]
            local θ = plan.weight_right[e]
            for b = 1:Nb, a = 1:Nb
                shifted[a, b] =
                    l == r || θ == 0 ? G[l, m, a, b] :
                    (1 - θ) * G[l, m, a, b] + θ * G[r, m, a, b]
            end
            mul!(temporary, T, shifted)
            mul!(result, temporary, adjoint(T))
            for b = 1:Nb, a = 1:Nb
                Σ[e, m, a, b] = result[a, b]
            end
        end
    end
    return Σ
end

function _embedding_self_energy_production(
    problem::NEGFProblem,
    green::GreenState,
    cache::ProductionCache,
    options::ProductionOptions = ProductionOptions(),
)
    bp = problem.basis
    plusR = _embedding_component_production(bp.T₊, green.Gᴿ, cache.W̃₊ᴱᵖ, options)
    plusL = _embedding_component_production(bp.T₊, green.Gˡ, cache.W̃₊ᴱᵖ, options)
    plusG = _embedding_component_production(bp.T₊, green.Gᵍ, cache.W̃₊ᴱᵖ, options)
    minusR = _embedding_component_production(bp.T₋, green.Gᴿ, cache.W̃₋ᴱᵖ, options)
    minusL = _embedding_component_production(bp.T₋, green.Gˡ, cache.W̃₋ᴱᵖ, options)
    minusG = _embedding_component_production(bp.T₋, green.Gᵍ, cache.W̃₋ᴱᵖ, options)
    plus = SelfEnergyFamily(plusR, plusL, plusG)
    minus = SelfEnergyFamily(minusR, minusL, minusG)
    total = SelfEnergyFamily(plusR + minusR, plusL + minusL, plusG + minusG)
    return total, plus, minus
end

function _embedding_self_energy_selected(
    problem::NEGFProblem,
    green::GreenState,
    cache::ProductionCache,
    options::ProductionOptions;
    scattering = nothing,
    h = nothing,
    cavity_plan = nothing,
)
    if options.algorithms.embedding === :finite_chain
        scattering === nothing && throw(
            ArgumentError("finite_chain embedding requires irreducible local scattering"),
        )
        h === nothing && throw(
            ArgumentError(
                "finite_chain embedding requires the current Hartree Hamiltonian",
            ),
        )
        local_scattering = _sum_selfenergies(scattering, size(green.Gᴿ))
        return _cavity_embedding_production(
            problem,
            local_scattering,
            h,
            options;
            plan = cavity_plan,
        )
    end
    if options.algorithms.energy_shift in (:dense, :conservative_pair)
        return embedding_self_energy(problem, green)
    end
    return _embedding_self_energy_production(problem, green, cache, options)
end

function _selected_energy_shift(
    problem::NEGFProblem,
    cache::ProductionCache,
    direction::Symbol,
    X::AbstractArray{<:Number,4},
    options::ProductionOptions,
)
    destination = similar(X, promote_type(Float64, eltype(X)))
    return _selected_energy_shift!(destination, problem, cache, direction, X, options)
end

"""Apply the selected exact shift into independently owned output storage."""
function _selected_energy_shift!(
    destination::AbstractArray{<:Number,4},
    problem::NEGFProblem,
    cache::ProductionCache,
    direction::Symbol,
    X::AbstractArray{<:Number,4},
    options::ProductionOptions,
)
    size(destination) == size(X) || throw(DimensionMismatch("energy-shift output differs"))
    Base.mightalias(destination, X) &&
        throw(ArgumentError("energy-shift input and output must not alias"))
    if options.algorithms.energy_shift in (:dense, :conservative_pair)
        matrix =
            direction === :plus_lo ? problem.W₊ᴸᴼ :
            direction === :minus_lo ? problem.W₋ᴸᴼ :
            throw(ArgumentError("unsupported energy-shift direction $direction"))
        matrix isa EnergyShiftPlan &&
            return apply_energy_shift!(destination, matrix, X; options)
        NE, Nk, Nb, Nb2 = size(X)
        Nb == Nb2 && size(matrix) == (NE, NE) ||
            throw(DimensionMismatch("energy-shift matrix differs"))
        # Preserve the same column-wise matrix-vector operation as the reference.
        # Dense/conservative matrices and compact operators share this mul! API.
        for m = 1:Nk, a = 1:Nb, b = 1:Nb
            mul!(view(destination, :, m, a, b), matrix, view(X, :, m, a, b))
        end
        return destination
    end
    plan =
        direction === :plus_lo ? cache.W̃₊ᴸᴼ :
        direction === :minus_lo ? cache.W̃₋ᴸᴼ :
        throw(ArgumentError("unsupported energy-shift direction $direction"))
    return apply_energy_shift!(destination, plan, X; options)
end

function _selected_static_contraction(
    problem::NEGFProblem,
    cache::ProductionCache,
    mechanism::Symbol,
    G::AbstractArray{<:Number,4},
    qᴷ::Real,
    options::ProductionOptions,
)
    if options.algorithms.contraction === :literal
        return _static_contraction(problem.kernels.K[mechanism], G, problem.grids.wᵏ, qᴷ)
    end
    return production_static_contraction(
        cache.kernels[mechanism],
        G,
        problem.grids.wᵏ,
        qᴷ;
        energy_chunk = options.energy_chunk,
        parallel_backend = options.parallel_backend,
        worker_count = options.worker_count,
    )
end

function _production_lo_contraction(
    op::AbstractProductionKernel,
    green::GreenState,
    problem::NEGFProblem,
    cache::ProductionCache,
    qᴷ::Real,
    options::ProductionOptions,
)
    if problem.models.lo_population !== :thermal
        emission_in = _selected_static_contraction(
            problem,
            cache,
            :LO,
            _selected_energy_shift(problem, cache, :plus_lo, green.Gˡ, options),
            qᴷ,
            options,
        )
        absorption_in = _selected_static_contraction(
            problem,
            cache,
            :LO,
            _selected_energy_shift(problem, cache, :minus_lo, green.Gˡ, options),
            qᴷ,
            options,
        )
        emission_out = _selected_static_contraction(
            problem,
            cache,
            :LO,
            _selected_energy_shift(problem, cache, :minus_lo, green.Gᵍ, options),
            qᴷ,
            options,
        )
        absorption_out = _selected_static_contraction(
            problem,
            cache,
            :LO,
            _selected_energy_shift(problem, cache, :plus_lo, green.Gᵍ, options),
            qᴷ,
            options,
        )
        population = _lo_population(problem, green, emission_out, absorption_out)
        return (population+1) .* emission_in .+ population .* absorption_in,
        (population+1) .* emission_out .+ population .* absorption_out
    end
    Nᴸᴼ = _thermal_lo_occupation(problem.physical)
    shape = size(green.Gˡ)
    if options.algorithms.contraction === :literal
        Gˡ₊ = _selected_energy_shift(problem, cache, :plus_lo, green.Gˡ, options)
        Gˡ₋ = _selected_energy_shift(problem, cache, :minus_lo, green.Gˡ, options)
        Gᵍ₊ = _selected_energy_shift(problem, cache, :plus_lo, green.Gᵍ, options)
        Gᵍ₋ = _selected_energy_shift(problem, cache, :minus_lo, green.Gᵍ, options)
        Σˡ =
            (Nᴸᴼ + 1) .*
            _selected_static_contraction(problem, cache, :LO, Gˡ₊, qᴷ, options) .+
            Nᴸᴼ .* _selected_static_contraction(problem, cache, :LO, Gˡ₋, qᴷ, options)
        Σᵍ =
            (Nᴸᴼ + 1) .*
            _selected_static_contraction(problem, cache, :LO, Gᵍ₋, qᴷ, options) .+
            Nᴸᴼ .* _selected_static_contraction(problem, cache, :LO, Gᵍ₊, qᴷ, options)
        return Σˡ, Σᵍ
    end
    Σˡ = zeros(ComplexF64, shape)
    shifted = _selected_energy_shift(problem, cache, :plus_lo, green.Gˡ, options)
    _production_static_contraction!(
        Σˡ,
        op,
        shifted,
        problem.grids.wᵏ,
        qᴷ,
        options.energy_chunk;
        α = Nᴸᴼ + 1,
        β = 0,
        parallel_backend = options.parallel_backend,
        worker_count = options.worker_count,
    )
    _selected_energy_shift!(shifted, problem, cache, :minus_lo, green.Gˡ, options)
    _production_static_contraction!(
        Σˡ,
        op,
        shifted,
        problem.grids.wᵏ,
        qᴷ,
        options.energy_chunk;
        α = Nᴸᴼ,
        β = 1,
        parallel_backend = options.parallel_backend,
        worker_count = options.worker_count,
    )

    Σᵍ = zeros(ComplexF64, shape)
    _selected_energy_shift!(shifted, problem, cache, :minus_lo, green.Gᵍ, options)
    _production_static_contraction!(
        Σᵍ,
        op,
        shifted,
        problem.grids.wᵏ,
        qᴷ,
        options.energy_chunk;
        α = Nᴸᴼ + 1,
        β = 0,
        parallel_backend = options.parallel_backend,
        worker_count = options.worker_count,
    )
    _selected_energy_shift!(shifted, problem, cache, :plus_lo, green.Gᵍ, options)
    _production_static_contraction!(
        Σᵍ,
        op,
        shifted,
        problem.grids.wᵏ,
        qᴷ,
        options.energy_chunk;
        α = Nᴸᴼ,
        β = 1,
        parallel_backend = options.parallel_backend,
        worker_count = options.worker_count,
    )
    return Σˡ, Σᵍ
end

function _retarded_self_energy_selected(
    Σˡ,
    Σᵍ,
    problem::NEGFProblem,
    cache::ProductionCache,
    options::ProductionOptions;
    return_roundoff::Bool = false,
)
    algorithms = options.algorithms
    if algorithms.retarded_real_part === :drop
        Γ = ComplexF64.(im .* (Σᵍ .- Σˡ))
        Σᴿ = ComplexF64.(-0.5im .* Γ)
        return return_roundoff ? (Σᴿ, Γ, 0.0) : (Σᴿ, Γ)
    elseif algorithms.hilbert === :product_integration
        Γ = ComplexF64.(im .* (Σᵍ .- Σˡ))
        Λ, residual = product_integration_hilbert_transform(
            Γ,
            problem.grids.ε,
            problem.grids.wᴱ;
            return_roundoff = true,
            operator = cache.product_hilbert_operator,
        )
        Σᴿ = ComplexF64.(Λ .- 0.5im .* Γ)
        return return_roundoff ? (Σᴿ, Γ, residual) : (Σᴿ, Γ)
    elseif algorithms.hilbert === :direct
        return retarded_self_energy(Σˡ, Σᵍ, problem.grids; return_roundoff)
    end
    return _retarded_self_energy_fft_production(
        Σˡ,
        Σᵍ,
        problem.grids,
        cache.fft_hilbert_plan,
        options;
        return_roundoff,
    )
end

function _diagonalize_blocks!(X::AbstractArray{ComplexF64,4})
    _, _, Nb, Nb2 = size(X)
    Nb == Nb2 || throw(DimensionMismatch("self-energy blocks must be square"))
    for b = 1:Nb, a = 1:Nb
        a == b && continue
        @views fill!(X[:, :, a, b], 0)
    end
    return X
end

function _average_transverse_momentum!(
    X::AbstractArray{ComplexF64,4},
    weights::AbstractVector{<:Real},
)
    NE, Nk, Nb, Nb2 = size(X)
    Nb == Nb2 || throw(DimensionMismatch("self-energy blocks must be square"))
    length(weights) == Nk || throw(DimensionMismatch("radial weights differ"))
    normalization = sum(weights)
    normalization > 0 || throw(ArgumentError("radial weights must have positive sum"))
    for b = 1:Nb, a = 1:Nb, e = 1:NE
        value = zero(ComplexF64)
        for m = 1:Nk
            value += weights[m] * X[e, m, a, b]
        end
        value /= normalization
        for m = 1:Nk
            X[e, m, a, b] = value
        end
    end
    return X
end

function _apply_scattering_model!(
    family::SelfEnergyFamily,
    problem::NEGFProblem,
    algorithms::AlgorithmOptions,
)
    if algorithms.self_energy_structure === :diagonal
        _diagonalize_blocks!(family.Σᴿ)
        _diagonalize_blocks!(family.Σˡ)
        _diagonalize_blocks!(family.Σᵍ)
    end
    if algorithms.transverse_momentum === :averaged
        _average_transverse_momentum!(family.Σᴿ, problem.grids.wᵏ)
        _average_transverse_momentum!(family.Σˡ, problem.grids.wᵏ)
        _average_transverse_momentum!(family.Σᵍ, problem.grids.wᵏ)
    end
    return family
end

function _scattering_candidate_production(
    problem::NEGFProblem,
    green::GreenState,
    cache::ProductionCache,
    options::ProductionOptions;
    return_roundoff::Bool = false,
)
    candidates = Dict{Symbol,SelfEnergyFamily}()
    roundoff = 0.0
    for mechanism in problem.kernels.enabled
        if mechanism === :electron_electron
            components = sppa_components(problem, green)
            sigmaR, _, residual = _retarded_self_energy_selected(
                components.lesser,
                components.greater,
                problem,
                cache,
                options;
                return_roundoff = true,
            )
            for m in axes(sigmaR, 2), e in axes(sigmaR, 1)
                @views sigmaR[e, m, :, :] .+= components.exchange
            end
            candidates[mechanism] = _apply_scattering_model!(
                SelfEnergyFamily(sigmaR, components.lesser, components.greater),
                problem,
                options.algorithms,
            )
            roundoff = max(roundoff, residual)
            continue
        end
        op = cache.kernels[mechanism]
        qᴷ = problem.kernels.qᴷ[mechanism]
        if mechanism === :LO
            Σˡ, Σᵍ = _production_lo_contraction(op, green, problem, cache, qᴷ, options)
        else
            Σˡ = _selected_static_contraction(
                problem,
                cache,
                mechanism,
                green.Gˡ,
                qᴷ,
                options,
            )
            Σᵍ = _selected_static_contraction(
                problem,
                cache,
                mechanism,
                green.Gᵍ,
                qᴷ,
                options,
            )
        end
        if return_roundoff && options.verify_fft_roundoff
            Σᴿ, _, residual = _retarded_self_energy_selected(
                Σˡ,
                Σᵍ,
                problem,
                cache,
                options;
                return_roundoff = true,
            )
            roundoff = max(roundoff, residual)
        else
            Σᴿ, _ = _retarded_self_energy_selected(Σˡ, Σᵍ, problem, cache, options)
        end
        family = SelfEnergyFamily(Σᴿ, Σˡ, Σᵍ)
        candidates[mechanism] =
            _apply_scattering_model!(family, problem, options.algorithms)
    end
    return return_roundoff ? (candidates, roundoff) : candidates
end

function _retarded_green_production(
    ε::AbstractVector{<:Real},
    h::AbstractArray{<:Number,3},
    Σᴿ::AbstractArray{<:Number,4};
    η::Real = 0.0,
    options::ProductionOptions = ProductionOptions(),
)
    NE, Nk, Nb, Nb2 = size(Σᴿ)
    Nb == Nb2 || throw(DimensionMismatch("Σᴿ blocks must be square"))
    size(h) == (Nk, Nb, Nb) || throw(DimensionMismatch("h has wrong shape"))
    length(ε) == NE || throw(DimensionMismatch("energy axis differs"))
    Gᴿ = zeros(ComplexF64, size(Σᴿ))
    κD = zeros(Float64, NE, Nk)
    scaleD = zeros(Float64, NE, Nk)
    blocks = NE * Nk
    workers = _production_worker_count(options, blocks)
    _production_worker_ranges!(blocks, workers) do _, first_linear, last_linear
        local D = zeros(ComplexF64, Nb, Nb)
        local inverse = similar(D)
        local inverse_workspace = _ProductionInverseWorkspace(Nb)
        @inbounds for linear = first_linear:last_linear
            local e = (linear - 1) % NE + 1
            local m = (linear - 1) ÷ NE + 1
            for b = 1:Nb, a = 1:Nb
                D[a, b] =
                    (a == b ? ε[e] + im * η : 0.0 + 0.0im) - h[m, a, b] - Σᴿ[e, m, a, b]
            end
            local scale = 0.0
            for a = 1:Nb
                local row_sum = 0.0
                for b = 1:Nb
                    row_sum += abs(D[a, b])
                end
                scale = max(scale, row_sum)
            end
            scale > 0 || throw(SingularException(0))
            for b = 1:Nb, a = 1:Nb
                inverse[a, b] = D[a, b] / scale
            end
            _production_inverse!(inverse, inverse_workspace)
            for b = 1:Nb, a = 1:Nb
                Gᴿ[e, m, a, b] = inverse[a, b] / scale
            end
            κD[e, m] = cond(D)
            scaleD[e, m] = scale
        end
    end
    return Gᴿ, κD, scaleD
end

function _green_state_production(
    Gᴿ::AbstractArray{<:Number,4},
    Gˡ::AbstractArray{<:Number,4},
    κD,
    scaleD,
    options::ProductionOptions = ProductionOptions();
    greater::Union{Nothing,Array{ComplexF64,4}} = nothing,
)
    size(Gᴿ) == size(Gˡ) || throw(DimensionMismatch("Green arrays differ"))
    NE, Nk = size(Gᴿ, 1), size(Gᴿ, 2)
    Nb = size(Gᴿ, 3)
    greater === nothing ||
        size(greater) == size(Gᴿ) ||
        throw(DimensionMismatch("greater Green array differs"))
    Gᵍ = greater === nothing ? zeros(ComplexF64, size(Gᴿ)) : greater
    A = zeros(ComplexF64, size(Gᴿ))
    _production_parallel_linear!(NE * Nk, options) do linear
        local e = (linear - 1) % NE + 1
        local m = (linear - 1) ÷ NE + 1
        @inbounds for b = 1:Nb, a = 1:Nb
            local difference = Gᴿ[e, m, a, b] - conj(Gᴿ[e, m, b, a])
            greater === nothing && (Gᵍ[e, m, a, b] = difference + Gˡ[e, m, a, b])
            A[e, m, a, b] = im * difference
        end
    end
    GR_array = Gᴿ isa Array{ComplexF64,4} ? Gᴿ : ComplexF64.(Gᴿ)
    GL_array = Gˡ isa Array{ComplexF64,4} ? Gˡ : ComplexF64.(Gˡ)
    return GreenState(GR_array, GL_array, Gᵍ, A, κD, scaleD)
end

function _number_functional_production(
    Gˡ::AbstractArray{<:Number,4},
    grids::ModelGrids,
    g_s::Integer,
    options::ProductionOptions,
)
    NE, Nk, Nb, Nb2 = size(Gˡ)
    Nb == Nb2 || throw(DimensionMismatch("Gˡ blocks must be square"))
    length(grids.wᴱ) == NE && length(grids.wᵏ) == Nk ||
        throw(DimensionMismatch("quadrature grids differ from Gˡ"))
    terms = Vector{ComplexF64}(undef, NE * Nk)
    _production_parallel_linear!(length(terms), options) do linear
        e = (linear - 1) ÷ Nk + 1
        m = (linear - 1) % Nk + 1
        block_trace = 0.0 + 0.0im
        @inbounds for a = 1:Nb
            block_trace += Gˡ[e, m, a, a]
        end
        terms[linear] = grids.wᴱ[e] * grids.wᵏ[m] * block_trace
    end
    value = 0.0 + 0.0im
    @inbounds for term in terms
        value += term
    end
    value *= -im * g_s / (2π)
    abs(imag(value)) ≤ _NUMBER_FUNCTIONAL_REALNESS_LIMIT * max(abs(real(value)), 1.0) ||
        throw(DomainError(value, "number functional is not real"))
    return real(value)
end

function _total_family_production(
    scattering::Dict{Symbol,SelfEnergyFamily},
    embedding::SelfEnergyFamily,
    mechanism_order,
    options::ProductionOptions,
)
    families = [scattering[name] for name in mechanism_order]
    total = SelfEnergyFamily(
        similar(embedding.Σᴿ),
        similar(embedding.Σˡ),
        similar(embedding.Σᵍ),
    )
    _production_parallel_linear!(length(embedding.Σᴿ), options) do index
        valueR = embedding.Σᴿ[index]
        valueL = embedding.Σˡ[index]
        valueG = embedding.Σᵍ[index]
        for family in families
            valueR += family.Σᴿ[index]
            valueL += family.Σˡ[index]
            valueG += family.Σᵍ[index]
        end
        total.Σᴿ[index] = valueR
        total.Σˡ[index] = valueL
        total.Σᵍ[index] = valueG
    end
    return total
end

function _mix_family_production!(
    old::SelfEnergyFamily,
    candidate::SelfEnergyFamily,
    α::Real,
    options::ProductionOptions,
)
    0 < α ≤ 1 || throw(ArgumentError("mixing factor must lie in (0,1]"))
    length(old.Σᴿ) == length(candidate.Σᴿ) ||
        throw(DimensionMismatch("self-energy families differ"))
    β = 1 - α
    _production_parallel_linear!(length(old.Σᴿ), options) do index
        old.Σᴿ[index] = β * old.Σᴿ[index] + α * candidate.Σᴿ[index]
        old.Σˡ[index] = β * old.Σˡ[index] + α * candidate.Σˡ[index]
        old.Σᵍ[index] = β * old.Σᵍ[index] + α * candidate.Σᵍ[index]
    end
    return old
end

function _mix_mechanisms_production!(
    old::Dict{Symbol,SelfEnergyFamily},
    candidate::Dict{Symbol,SelfEnergyFamily},
    mechanism_order,
    α::Real,
    options::ProductionOptions,
)
    for name in mechanism_order
        _mix_family_production!(old[name], candidate[name], α, options)
    end
    return old
end

function _embedding_sum_production!(
    embedding::SelfEnergyFamily,
    plus::SelfEnergyFamily,
    minus::SelfEnergyFamily,
    options::ProductionOptions,
)
    _production_parallel_linear!(length(embedding.Σᴿ), options) do index
        embedding.Σᴿ[index] = plus.Σᴿ[index] + minus.Σᴿ[index]
        embedding.Σˡ[index] = plus.Σˡ[index] + minus.Σˡ[index]
        embedding.Σᵍ[index] = plus.Σᵍ[index] + minus.Σᵍ[index]
    end
    return embedding
end
