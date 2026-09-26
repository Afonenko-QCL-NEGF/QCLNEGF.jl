export ProductionResidualWorkspace,
    production_residual_suite, production_selfenergy_residual

"""
Reusable storage for the exact production SCBA residuals.

The four per-block vectors make the two global Frobenius-norm reductions
independent of thread scheduling: workers write one scalar per `(E,k)` block,
and the final sum is performed serially in canonical energy-major order.  The
small matrices are private to a Julia worker.  Consequently, the threaded
implementation changes neither the residual definitions nor their reduction
order.

See [the SCBA fixed point](@ref theory-scba),
[verification residuals](@ref theory-validation), and
[the production backend](@ref theory-production).
"""
struct ProductionResidualWorkspace
    N_E::Int
    N_k::Int
    N_b::Int
    chunk_size::Int
    spectral_numerator::Vector{Float64}
    spectral_denominator::Vector{Float64}
    keldysh_numerator::Vector{Float64}
    keldysh_denominator::Vector{Float64}
    dyson_partial::Vector{Float64}
    positivity_partial::Vector{Float64}
    positivity_witness::Vector{_PositivityWitness}
    causality_partial::Vector{Float64}
    D::Vector{Matrix{ComplexF64}}
    Γ::Vector{Matrix{ComplexF64}}
    temporary_1::Vector{Matrix{ComplexF64}}
    temporary_2::Vector{Matrix{ComplexF64}}
    hermitian::Vector{Matrix{ComplexF64}}
    lapack_values::Vector{Vector{Float64}}
    lapack_work::Vector{Vector{ComplexF64}}
    lapack_rwork::Vector{Vector{Float64}}
    lapack_info::Vector{Vector{LinearAlgebra.BlasInt}}
end

"""
    ProductionResidualWorkspace(N_E, N_k, N_b;
                                chunk_size=64, worker_count=0)

Allocate the bounded workspace used by [`production_residual_suite`](@ref)
and [`production_selfenergy_residual`](@ref).  `chunk_size` controls only task
granularity.  It does not change the canonical order used for the two sums.
`worker_count=0` selects the complete default Julia thread pool; use
`worker_count=1` with a BLAS-parallel production backend.  One workspace is
owned by one solver call and must not be shared by concurrent calls.
"""
function ProductionResidualWorkspace(
    N_E::Integer,
    N_k::Integer,
    N_b::Integer;
    chunk_size::Integer = 64,
    worker_count::Integer = 0,
)
    N_E > 0 || throw(ArgumentError("N_E must be positive"))
    N_k > 0 || throw(ArgumentError("N_k must be positive"))
    N_b > 0 || throw(ArgumentError("N_b must be positive"))
    chunk_size > 0 || throw(ArgumentError("chunk_size must be positive"))
    worker_count ≥ 0 || throw(ArgumentError("worker_count cannot be negative"))
    NE, Nk, Nb = Int(N_E), Int(N_k), Int(N_b)
    blocks = Base.checked_mul(NE, Nk)
    chunks = cld(blocks, Int(chunk_size))
    # Logical workers, not absolute thread IDs.  With `--threads N,1`, default
    # pool IDs are offset by the interactive pool; indexing by `threadid()` is
    # therefore incorrect.  Explicit spawned tasks below own slots `1:workers`.
    available = Base.Threads.nthreads(:default)
    workers = worker_count == 0 ? available : min(Int(worker_count), available)
    scalar_blocks() = zeros(Float64, blocks)
    scalar_chunks() = zeros(Float64, chunks)
    matrices() = [zeros(ComplexF64, Nb, Nb) for _ = 1:workers]
    real_vectors(length) = [zeros(Float64, length) for _ = 1:workers]
    complex_vectors(length) = [zeros(ComplexF64, length) for _ = 1:workers]
    integer_vectors(length) = [zeros(LinearAlgebra.BlasInt, length) for _ = 1:workers]
    return ProductionResidualWorkspace(
        NE,
        Nk,
        Nb,
        Int(chunk_size),
        scalar_blocks(),
        scalar_blocks(),
        scalar_blocks(),
        scalar_blocks(),
        scalar_chunks(),
        scalar_chunks(),
        fill(_empty_positivity_witness(), chunks),
        scalar_chunks(),
        matrices(),
        matrices(),
        matrices(),
        matrices(),
        matrices(),
        real_vectors(Nb),
        complex_vectors(max(1, 3Nb)),
        real_vectors(max(1, 5Nb)),
        integer_vectors(1),
    )
end

ProductionResidualWorkspace(
    problem::NEGFProblem;
    chunk_size::Integer = 64,
    worker_count::Integer = 0,
) = ProductionResidualWorkspace(
    problem.numerical.N_E,
    problem.numerical.N_k,
    problem.numerical.N_b;
    chunk_size,
    worker_count,
)

@inline _production_residual_chunks(workspace::ProductionResidualWorkspace) =
    length(workspace.dyson_partial)

# Allocated capacity and active phase width are separate. A smaller resource
# lease reuses the same private buffers without changing canonical reductions.
@inline function _production_residual_worker_count(workspace, chunks, requested::Integer)
    requested >= 0 || throw(ArgumentError("worker_count cannot be negative"))
    capacity = min(length(workspace.D), chunks, Base.Threads.nthreads(:default))
    return requested == 0 ? capacity : min(Int(requested), capacity)
end

@inline function _production_block_indices(linear::Int, N_k::Int)
    # Match validation.jl exactly: `for e in 1:N_E, m in 1:N_k`.
    e = (linear - 1) ÷ N_k + 1
    m = (linear - 1) % N_k + 1
    return e, m
end

@inline function _production_frobenius_squared(X::AbstractMatrix)
    value = 0.0
    @inbounds for j in axes(X, 2), i in axes(X, 1)
        value += abs2(X[i, j])
    end
    return value
end

@inline function _production_frobenius_norm(X::AbstractMatrix)
    scale = 0.0
    sum_squares = 1.0
    @inbounds for j in axes(X, 2), i in axes(X, 1)
        magnitude = abs(X[i, j])
        isfinite(magnitude) || return Inf
        if magnitude != 0
            if scale < magnitude
                sum_squares = 1 + sum_squares * (scale / magnitude)^2
                scale = magnitude
            else
                sum_squares += (magnitude / scale)^2
            end
        end
    end
    return scale == 0 ? 0.0 : scale * sqrt(sum_squares)
end

# The `(E,k,:,: )` blocks are strided, non-contiguous views.  Sending them to
# the generic `mul!` path makes Julia materialize BLAS temporaries for every
# five-by-five product.  These canonical small-block loops are allocation-free
# and evaluate the same dense products directly from the stored arrays.
@inline function _production_matrix_product!(
    C::Matrix{ComplexF64},
    A::AbstractMatrix,
    B::AbstractMatrix,
)
    Nb = size(C, 1)
    @inbounds for j = 1:Nb, i = 1:Nb
        value = 0.0 + 0.0im
        for k = 1:Nb
            value += A[i, k] * B[k, j]
        end
        C[i, j] = value
    end
    return C
end

@inline function _production_matrix_product_adjoint_right!(
    C::Matrix{ComplexF64},
    A::AbstractMatrix,
    B::AbstractMatrix,
)
    Nb = size(C, 1)
    @inbounds for j = 1:Nb, i = 1:Nb
        value = 0.0 + 0.0im
        for k = 1:Nb
            value += A[i, k] * conj(B[j, k])
        end
        C[i, j] = value
    end
    return C
end

@inline function _production_allfinite(X::AbstractMatrix)
    @inbounds for j in axes(X, 2), i in axes(X, 1)
        isfinite(X[i, j]) || return false
    end
    return true
end

@inline _production_metric(value::Real) = isfinite(value) ? Float64(value) : Inf

function _check_production_residual_workspace(
    workspace::ProductionResidualWorkspace,
    N_E::Int,
    N_k::Int,
    N_b::Int,
)
    (workspace.N_E, workspace.N_k, workspace.N_b) == (N_E, N_k, N_b) ||
        throw(DimensionMismatch("production residual workspace has wrong dimensions"))
    workspace.chunk_size > 0 ||
        throw(DimensionMismatch("production residual workspace has invalid chunks"))
    blocks = Base.checked_mul(N_E, N_k)
    chunks = cld(blocks, workspace.chunk_size)
    for values in (
        workspace.spectral_numerator,
        workspace.spectral_denominator,
        workspace.keldysh_numerator,
        workspace.keldysh_denominator,
    )
        length(values) == blocks || throw(
            DimensionMismatch("production residual workspace has invalid block storage"),
        )
    end
    for values in (
        workspace.dyson_partial,
        workspace.positivity_partial,
        workspace.positivity_witness,
        workspace.causality_partial,
    )
        length(values) == chunks || throw(
            DimensionMismatch("production residual workspace has invalid chunk storage"),
        )
    end
    workers = length(workspace.D)
    workers > 0 ||
        throw(DimensionMismatch("production residual workspace has no worker buffers"))
    worker_fields = (
        workspace.Γ,
        workspace.temporary_1,
        workspace.temporary_2,
        workspace.hermitian,
        workspace.lapack_values,
        workspace.lapack_work,
        workspace.lapack_rwork,
        workspace.lapack_info,
    )
    all(length(field) == workers for field in worker_fields) ||
        throw(DimensionMismatch("production residual workspace worker counts differ"))
    for worker = 1:workers
        for matrix in (
            workspace.D[worker],
            workspace.Γ[worker],
            workspace.temporary_1[worker],
            workspace.temporary_2[worker],
            workspace.hermitian[worker],
        )
            size(matrix) == (N_b, N_b) || throw(
                DimensionMismatch("production residual workspace matrix has wrong size"),
            )
        end
        length(workspace.lapack_values[worker]) ≥ N_b ||
            throw(DimensionMismatch("production LAPACK values buffer is too short"))
        length(workspace.lapack_work[worker]) ≥ max(1, 3N_b) ||
            throw(DimensionMismatch("production LAPACK work buffer is too short"))
        length(workspace.lapack_rwork[worker]) ≥ max(1, 5N_b) ||
            throw(DimensionMismatch("production LAPACK rwork buffer is too short"))
        !isempty(workspace.lapack_info[worker]) ||
            throw(DimensionMismatch("production LAPACK info buffer is empty"))
    end
    return workspace
end

@inline function _resolve_production_residual_workspace(
    workspace,
    N_E::Int,
    N_k::Int,
    N_b::Int,
    chunk_size::Integer,
)
    if workspace === nothing
        return ProductionResidualWorkspace(N_E, N_k, N_b; chunk_size)
    end
    return _check_production_residual_workspace(workspace, N_E, N_k, N_b)
end

function _check_production_residual_inputs(
    ε,
    h,
    green::GreenState,
    total::SelfEnergyFamily,
    candidate::SelfEnergyFamily,
)
    NE, Nk, Nb, Nb2 = size(green.Gᴿ)
    Nb == Nb2 || throw(DimensionMismatch("Green blocks must be square"))
    for name in (:Gˡ, :Gᵍ, :A)
        size(getfield(green, name)) == (NE, Nk, Nb, Nb) ||
            throw(DimensionMismatch("Green arrays have different shapes"))
    end
    length(ε) == NE || throw(DimensionMismatch("energy axis differs"))
    size(h) == (Nk, Nb, Nb) || throw(DimensionMismatch("h has wrong shape"))
    for family in (total, candidate), name in (:Σᴿ, :Σˡ, :Σᵍ)
        size(getfield(family, name)) == (NE, Nk, Nb, Nb) ||
            throw(DimensionMismatch("self-energy array has wrong shape"))
    end
    return NE, Nk, Nb
end

"""Allocation-free LAPACK `zheev` using storage owned by one worker."""
function _production_hermitian_eigvals!(
    H::Matrix{ComplexF64},
    values::Vector{Float64},
    work::Vector{ComplexF64},
    rwork::Vector{Float64},
    info::Vector{LinearAlgebra.BlasInt},
)
    n = LinearAlgebra.BlasInt(size(H, 1))
    lwork = LinearAlgebra.BlasInt(length(work))
    info[1] = 0
    ccall(
        (LinearAlgebra.BLAS.@blasfunc(zheev_), LinearAlgebra.LAPACK.libblastrampoline),
        Cvoid,
        (
            Ref{UInt8},
            Ref{UInt8},
            Ref{LinearAlgebra.BlasInt},
            Ptr{ComplexF64},
            Ref{LinearAlgebra.BlasInt},
            Ptr{Float64},
            Ptr{ComplexF64},
            Ref{LinearAlgebra.BlasInt},
            Ptr{Float64},
            Ptr{LinearAlgebra.BlasInt},
            Clong,
            Clong,
        ),
        'N',
        'U',
        n,
        H,
        n,
        values,
        work,
        lwork,
        rwork,
        info,
        1,
        1,
    )
    return info[1] == 0
end

"""Allocation-free LAPACK `zgesvd` spectral norm using worker storage."""
function _production_opnorm!(
    A::Matrix{ComplexF64},
    values::Vector{Float64},
    work::Vector{ComplexF64},
    rwork::Vector{Float64},
    info::Vector{LinearAlgebra.BlasInt},
)
    n = LinearAlgebra.BlasInt(size(A, 1))
    lwork = LinearAlgebra.BlasInt(length(work))
    one = LinearAlgebra.BlasInt(1)
    info[1] = 0
    # With JOBU=JOBVT='N', LAPACK does not reference U or VT.  `work` is a
    # valid non-null dummy pointer for both formal arguments.
    ccall(
        (LinearAlgebra.BLAS.@blasfunc(zgesvd_), LinearAlgebra.LAPACK.libblastrampoline),
        Cvoid,
        (
            Ref{UInt8},
            Ref{UInt8},
            Ref{LinearAlgebra.BlasInt},
            Ref{LinearAlgebra.BlasInt},
            Ptr{ComplexF64},
            Ref{LinearAlgebra.BlasInt},
            Ptr{Float64},
            Ptr{ComplexF64},
            Ref{LinearAlgebra.BlasInt},
            Ptr{ComplexF64},
            Ref{LinearAlgebra.BlasInt},
            Ptr{ComplexF64},
            Ref{LinearAlgebra.BlasInt},
            Ptr{Float64},
            Ptr{LinearAlgebra.BlasInt},
            Clong,
            Clong,
        ),
        'N',
        'N',
        n,
        n,
        A,
        n,
        values,
        work,
        one,
        work,
        one,
        work,
        lwork,
        rwork,
        info,
        1,
        1,
    )
    info[1] == 0 || return Inf
    return _production_metric(maximum(values))
end

"""Fill `H` with the Hermitian part of `factor*X` and return its PSD defect."""
function _production_negative_part!(
    H::Matrix{ComplexF64},
    X::AbstractMatrix,
    factor::ComplexF64,
    values::Vector{Float64},
    work::Vector{ComplexF64},
    rwork::Vector{Float64},
    info::Vector{LinearAlgebra.BlasInt};
    floor::Float64 = 1e-14,
    witness_sink = nothing,
    witness_chunk::Int = 0,
    matrix_kind::Symbol = :none,
    energy_index::Int = 0,
    momentum_index::Int = 0,
)
    _production_allfinite(X) || return Inf
    Nb = size(H, 1)
    @inbounds for j = 1:Nb, i = 1:Nb
        H[i, j] = (factor * X[i, j] + conj(factor * X[j, i])) / 2
    end
    _production_allfinite(H) || return Inf
    _production_hermitian_eigvals!(H, values, work, rwork, info) || return Inf
    all(isfinite, values) || return Inf
    hermiticity_squared = 0.0
    @inbounds for j = 1:Nb, i = 1:Nb
        hermiticity_squared += abs2(factor*X[i, j] - conj(factor*X[j, i]))
    end
    witness = _positivity_witness(
        matrix_kind,
        energy_index,
        momentum_index,
        values;
        floor,
        hermiticity_defect = sqrt(hermiticity_squared)/2,
    )
    if witness_sink !== nothing && witness.ratio > witness_sink[witness_chunk].ratio
        witness_sink[witness_chunk] = witness
    end
    return _production_metric(witness.ratio)
end

"""Fill `H` with `Im(X)` in the retarded-matrix sense."""
function _production_retarded_imaginary!(H::Matrix{ComplexF64}, X::AbstractMatrix)
    Nb = size(H, 1)
    @inbounds for j = 1:Nb, i = 1:Nb
        H[i, j] = (X[i, j] - conj(X[j, i])) / (2im)
    end
    return H
end

function _production_causality_block!(
    H::Matrix{ComplexF64},
    norm_scratch::Matrix{ComplexF64},
    ΣR::AbstractMatrix,
    values::Vector{Float64},
    work::Vector{ComplexF64},
    rwork::Vector{Float64},
    info::Vector{LinearAlgebra.BlasInt};
    floor::Float64 = 1e-14,
)
    _production_allfinite(ΣR) || return Inf
    _production_retarded_imaginary!(H, ΣR)
    _production_allfinite(H) || return Inf
    _production_hermitian_eigvals!(H, values, work, rwork, info) || return Inf
    all(isfinite, values) || return Inf
    numerator = max(0.0, maximum(values))
    # A finite Frobenius norm bounds the spectral norm. For an already causal
    # block the numerator is exactly zero, so its expensive SVD cannot change
    # the residual. If the bound overflows, retain the reference SVD policy.
    if iszero(numerator) && isfinite(_production_frobenius_norm(ΣR))
        return 0.0
    end
    # The reference definition normalizes by the spectral norm of the full,
    # generally non-normal retarded self-energy rather than by Im(Σᴿ).
    copyto!(norm_scratch, ΣR)
    normΣ = _production_opnorm!(norm_scratch, values, work, rwork, info)
    isfinite(normΣ) || return Inf
    denominator = normΣ + floor
    return _production_metric(numerator / denominator)
end

@inline function _canonical_sum(values::Vector{Float64})
    result = 0.0
    @inbounds for i in eachindex(values)
        value = values[i]
        isfinite(value) || return Inf
        result += value
        isfinite(result) || return Inf
    end
    return result
end

@inline function _canonical_max(values::Vector{Float64})
    result = 0.0
    @inbounds for i in eachindex(values)
        value = values[i]
        isfinite(value) || return Inf
        result = max(result, value)
    end
    return result
end

"""
    production_residual_suite(ε, h, green, total, candidate;
                              workspace=nothing, chunk_size=64, worker_count=0)

Evaluate the exact production residual tuple
`(r_D, r_A, r_K, r_PSD, r_caus)` in one threaded pass over independent
`(E_e,k_m)` blocks.  These are algebraically identical to
`_dyson_residual`, `_spectral_residual`, `_keldysh_residual`,
`_green_positivity`, and `_causality_residual`.

Chunks have fixed boundaries and write to disjoint partials.  Frobenius terms
are stored per block and summed serially in energy-major order, making their
result independent of Julia thread scheduling and of `chunk_size`.  Any NaN,
infinity, eigensolver failure, or arithmetic overflow relevant to a metric is
reported as `Inf`, so a non-finite state can never satisfy a convergence test.

`worker_count` caps the active workers independently of allocated workspace
capacity; zero uses that capacity.

No approximation of the Dyson, Keldysh, spectral, positivity, or causality
equations is introduced.  See [verification residuals](@ref theory-validation)
and [the production SCBA loop](@ref theory-production).
"""
function production_residual_suite(
    ε::AbstractVector{<:Real},
    h::AbstractArray{<:Number,3},
    green::GreenState,
    total::SelfEnergyFamily,
    candidate::SelfEnergyFamily;
    workspace::Union{Nothing,ProductionResidualWorkspace} = nothing,
    chunk_size::Integer = 64,
    worker_count::Integer = 0,
)
    NE, Nk, Nb = _check_production_residual_inputs(ε, h, green, total, candidate)
    ws::ProductionResidualWorkspace =
        _resolve_production_residual_workspace(workspace, NE, Nk, Nb, chunk_size)
    blocks = NE * Nk
    chunks = _production_residual_chunks(ws)
    chunk = ws.chunk_size
    fill!(ws.dyson_partial, 0.0)
    fill!(ws.positivity_partial, 0.0)
    fill!(ws.positivity_witness, _empty_positivity_witness())
    fill!(ws.causality_partial, 0.0)

    workers = _production_residual_worker_count(ws, chunks, worker_count)
    @sync for worker = 1:workers
        Base.Threads.@spawn begin
            local D = ws.D[worker]
            local Γ = ws.Γ[worker]
            local temporary_1 = ws.temporary_1[worker]
            local temporary_2 = ws.temporary_2[worker]
            local H = ws.hermitian[worker]
            local values = ws.lapack_values[worker]
            local work = ws.lapack_work[worker]
            local rwork = ws.lapack_rwork[worker]
            local info = ws.lapack_info[worker]
            for c = worker:workers:chunks
                local rD = 0.0
                local rPSD = 0.0
                local rcaus = 0.0
                local first_linear = (c - 1) * chunk + 1
                local last_linear = min(c * chunk, blocks)
                @inbounds for linear = first_linear:last_linear
                    e, m = _production_block_indices(linear, Nk)
                    GR = view(green.Gᴿ, e, m, :, :)
                    GL = view(green.Gˡ, e, m, :, :)
                    GG = view(green.Gᵍ, e, m, :, :)
                    A = view(green.A, e, m, :, :)
                    ΣR = view(total.Σᴿ, e, m, :, :)
                    ΣL = view(total.Σˡ, e, m, :, :)
                    ΣG = view(total.Σᵍ, e, m, :, :)
                    candidate_ΣL = view(candidate.Σˡ, e, m, :, :)

                    # r_D = ||DᴿGᴿ-I||_F / (||Dᴿ||_F||Gᴿ||_F+sqrt(N_b)).
                    dyson_finite =
                        isfinite(ε[e]) &&
                        _production_allfinite(view(h, m, :, :)) &&
                        _production_allfinite(ΣR) &&
                        _production_allfinite(GR)
                    if dyson_finite
                        for j = 1:Nb, i = 1:Nb
                            D[i, j] = (i == j ? ε[e] : 0.0) - h[m, i, j] - ΣR[i, j]
                        end
                        _production_matrix_product!(temporary_1, D, GR)
                        for i = 1:Nb
                            temporary_1[i, i] -= 1
                        end
                        numerator = _production_frobenius_norm(temporary_1)
                        denominator =
                            _production_frobenius_norm(D) * _production_frobenius_norm(GR) +
                            sqrt(Nb)
                        value =
                            isfinite(numerator) && isfinite(denominator) ?
                            _production_metric(numerator / denominator) : Inf
                        rD = max(rD, value)
                    else
                        rD = Inf
                    end

                    # Γ=i(Σᵍ-Σˡ), used by both r_A and r_PSD without allocating a
                    # complete four-dimensional Γ array.
                    gamma_finite = _production_allfinite(ΣG) && _production_allfinite(ΣL)
                    if gamma_finite
                        for j = 1:Nb, i = 1:Nb
                            Γ[i, j] = im * (ΣG[i, j] - ΣL[i, j])
                        end
                    end

                    # Global spectral residual.  Each worker writes one pair of
                    # scalars; the reduction happens below in canonical order.
                    spectral_finite =
                        gamma_finite &&
                        _production_allfinite(GR) &&
                        _production_allfinite(A)
                    if spectral_finite
                        _production_matrix_product!(temporary_1, GR, Γ)
                        _production_matrix_product_adjoint_right!(
                            temporary_2,
                            temporary_1,
                            GR,
                        )
                        for j = 1:Nb, i = 1:Nb
                            temporary_2[i, j] = A[i, j] - temporary_2[i, j]
                        end
                        ws.spectral_numerator[linear] =
                            _production_metric(_production_frobenius_squared(temporary_2))
                        ws.spectral_denominator[linear] =
                            _production_metric(_production_frobenius_squared(A))
                    else
                        ws.spectral_numerator[linear] = Inf
                        ws.spectral_denominator[linear] = Inf
                    end

                    # Keldysh residual against the unmixed candidate self-energy.
                    keldysh_finite =
                        _production_allfinite(GR) &&
                        _production_allfinite(GL) &&
                        _production_allfinite(candidate_ΣL)
                    if keldysh_finite
                        _production_matrix_product!(temporary_1, GR, candidate_ΣL)
                        _production_matrix_product_adjoint_right!(
                            temporary_2,
                            temporary_1,
                            GR,
                        )
                        for j = 1:Nb, i = 1:Nb
                            temporary_2[i, j] = GL[i, j] - temporary_2[i, j]
                        end
                        ws.keldysh_numerator[linear] =
                            _production_metric(_production_frobenius_squared(temporary_2))
                        ws.keldysh_denominator[linear] =
                            _production_metric(_production_frobenius_squared(GL))
                    else
                        ws.keldysh_numerator[linear] = Inf
                        ws.keldysh_denominator[linear] = Inf
                    end

                    # Four positive-semidefinite matrices required by the NEGF sign
                    # convention.  Each eigensolve mutates only this worker's H.
                    for (kind, block, factor) in (
                        (:spectral, A, 1.0 + 0.0im),
                        (:occupied, GL, -1.0im),
                        (:unoccupied, GG, 1.0im),
                        (:broadening, Γ, 1.0 + 0.0im),
                    )
                        defect =
                            kind === :broadening && !gamma_finite ? Inf :
                            _production_negative_part!(
                                H,
                                block,
                                factor,
                                values,
                                work,
                                rwork,
                                info;
                                witness_sink = ws.positivity_witness,
                                witness_chunk = c,
                                matrix_kind = kind,
                                energy_index = e,
                                momentum_index = m,
                            )
                        if !isfinite(defect) && isfinite(ws.positivity_witness[c].ratio)
                            ws.positivity_witness[c] =
                                _PositivityWitness(kind, e, m, NaN, Inf, 1e-14, Inf)
                        end
                        rPSD = max(rPSD, defect)
                    end

                    # The fixed-number map must not hide an inadmissible raw
                    # Keldysh product. Test before redistribution as well as after.
                    for (kind, sigma, factor) in
                        ((:raw_occupied, ΣL, -1.0im), (:raw_unoccupied, ΣG, 1.0im))
                        _production_matrix_product!(temporary_1, GR, sigma)
                        _production_matrix_product_adjoint_right!(
                            temporary_2,
                            temporary_1,
                            GR,
                        )
                        defect = _production_negative_part!(
                            H,
                            temporary_2,
                            factor,
                            values,
                            work,
                            rwork,
                            info;
                            witness_sink = ws.positivity_witness,
                            witness_chunk = c,
                            matrix_kind = kind,
                            energy_index = e,
                            momentum_index = m,
                        )
                        rPSD = max(rPSD, defect)
                    end

                    # Retarded self-energy causality: λ_max(Im Σᴿ) <= 0.
                    rcaus = max(
                        rcaus,
                        _production_causality_block!(
                            H,
                            temporary_1,
                            ΣR,
                            values,
                            work,
                            rwork,
                            info,
                        ),
                    )
                end
                ws.dyson_partial[c] = rD
                ws.positivity_partial[c] = rPSD
                ws.causality_partial[c] = rcaus
            end
        end
    end

    spectral_numerator = _canonical_sum(ws.spectral_numerator)
    spectral_denominator = _canonical_sum(ws.spectral_denominator)
    keldysh_numerator = _canonical_sum(ws.keldysh_numerator)
    keldysh_denominator = _canonical_sum(ws.keldysh_denominator)
    rA =
        isfinite(spectral_numerator) && isfinite(spectral_denominator) ?
        _production_metric(
            sqrt(spectral_numerator) / (sqrt(spectral_denominator) + 1e-14),
        ) : Inf
    rK =
        isfinite(keldysh_numerator) && isfinite(keldysh_denominator) ?
        _production_metric(sqrt(keldysh_numerator) / (sqrt(keldysh_denominator) + 1e-14)) :
        Inf
    return (
        r_D = _canonical_max(ws.dyson_partial),
        r_A = rA,
        r_K = rK,
        r_PSD = _canonical_max(ws.positivity_partial),
        r_caus = _canonical_max(ws.causality_partial),
    )
end

"""Worst PSD block from the latest residual pass, without another eigenscan."""
function _production_positivity_witness(workspace::ProductionResidualWorkspace)
    result = _empty_positivity_witness()
    for witness in workspace.positivity_witness
        witness.ratio > result.ratio && (result = witness)
    end
    return result
end

@doc raw"""
    production_selfenergy_residual(old, candidate;
                                   workspace=nothing, chunk_size=64, worker_count=0)

Evaluate

```math
r_\Sigma=\max_{X\in\{R,<,>\}}
\frac{\|\Sigma^X_{\mathrm{candidate}}-\Sigma^X_{\mathrm{old}}\|_F}
     {\|\Sigma^X_{\mathrm{candidate}}\|_F+10^{-14}}
```

without allocating the four-dimensional difference arrays created by
`norm(candidate-old)`.  With a supplied `ProductionResidualWorkspace`, all
size-dependent storage is reused.  Per-block squared norms are written in
parallel and reduced in canonical order.  Non-finite data or overflow returns
`Inf`.

See [the SCBA convergence criteria](@ref theory-validation) and
[the production backend](@ref theory-production).
"""
function production_selfenergy_residual(
    old::SelfEnergyFamily,
    candidate::SelfEnergyFamily;
    workspace::Union{Nothing,ProductionResidualWorkspace} = nothing,
    chunk_size::Integer = 64,
    worker_count::Integer = 0,
)
    shape = size(candidate.Σᴿ)
    length(shape) == 4 || throw(DimensionMismatch("self-energy must be four-dimensional"))
    NE, Nk, Nb, Nb2 = shape
    Nb == Nb2 || throw(DimensionMismatch("self-energy blocks must be square"))
    for family in (old, candidate), name in (:Σᴿ, :Σˡ, :Σᵍ)
        size(getfield(family, name)) == shape ||
            throw(DimensionMismatch("self-energy arrays have different shapes"))
    end
    ws::ProductionResidualWorkspace =
        _resolve_production_residual_workspace(workspace, NE, Nk, Nb, chunk_size)
    blocks = NE * Nk
    chunks = _production_residual_chunks(ws)
    chunk = ws.chunk_size
    workers = _production_residual_worker_count(ws, chunks, worker_count)
    result = 0.0
    for name in (:Σᴿ, :Σˡ, :Σᵍ)
        Xold = getfield(old, name)
        Xnew = getfield(candidate, name)
        @sync for worker = 1:workers
            Base.Threads.@spawn begin
                for c = worker:workers:chunks
                    first_linear = (c - 1) * chunk + 1
                    last_linear = min(c * chunk, blocks)
                    @inbounds for linear = first_linear:last_linear
                        e, m = _production_block_indices(linear, Nk)
                        local difference_squared = 0.0
                        local candidate_squared = 0.0
                        local finite = true
                        for b = 1:Nb, a = 1:Nb
                            old_value = Xold[e, m, a, b]
                            new_value = Xnew[e, m, a, b]
                            if !(isfinite(old_value) && isfinite(new_value))
                                finite = false
                                break
                            end
                            difference_squared += abs2(new_value - old_value)
                            candidate_squared += abs2(new_value)
                        end
                        ws.spectral_numerator[linear] =
                            finite ? _production_metric(difference_squared) : Inf
                        ws.spectral_denominator[linear] =
                            finite ? _production_metric(candidate_squared) : Inf
                    end
                end
            end
        end
        difference_squared = _canonical_sum(ws.spectral_numerator)
        candidate_squared = _canonical_sum(ws.spectral_denominator)
        residual =
            isfinite(difference_squared) && isfinite(candidate_squared) ?
            _production_metric(
                sqrt(difference_squared) / (sqrt(candidate_squared) + 1e-14),
            ) : Inf
        result = max(result, residual)
    end
    return result
end

# Internal aliases keep the naming parallel to validation.jl and make the
# production call site explicit without adding private helpers to the public
# API documentation inventory.  The five projections deliberately share the
# same full signature: normal production code evaluates the suite once and
# destructures it, rather than traversing the arrays five times.
_production_residual_workspace(args...; kwargs...) =
    ProductionResidualWorkspace(args...; kwargs...)
_production_residual_suite(args...; kwargs...) =
    production_residual_suite(args...; kwargs...)
_dyson_residual_production(args...; kwargs...) =
    production_residual_suite(args...; kwargs...).r_D
_spectral_residual_production(args...; kwargs...) =
    production_residual_suite(args...; kwargs...).r_A
_keldysh_residual_production(args...; kwargs...) =
    production_residual_suite(args...; kwargs...).r_K
_green_positivity_production(args...; kwargs...) =
    production_residual_suite(args...; kwargs...).r_PSD
_causality_residual_production(args...; kwargs...) =
    production_residual_suite(args...; kwargs...).r_caus
_selfenergy_residual_production(args...; kwargs...) =
    production_selfenergy_residual(args...; kwargs...)
