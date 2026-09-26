"""
    build_shift_matrix(ε, δ)

Dense non-periodic piecewise-linear operator representing `X(ε+δ)`.
Rows whose shifted argument lies outside the energy window are exactly zero;
there is no wrap-around.  Shape: `(N_E,N_E)`.

Implements [EQ-SHIFT-001](@ref eq-energy-shift).
See [Field-periodic closure](@ref theory-periodicity) for the sign convention
and the reason the interpolation is non-cyclic.
"""
function build_shift_matrix(ε::AbstractVector{<:Real}, δ::Real)
    length(ε) >= 2 || throw(ArgumentError("at least two energy nodes are required"))
    all(isfinite, ε) && all(>(0), diff(ε)) ||
        throw(ArgumentError("energy nodes must be finite and strictly increasing"))
    isfinite(δ) || throw(ArgumentError("energy displacement must be finite"))
    N = length(ε)
    W = zeros(Float64, N, N)
    for e = 1:N
        target = ε[e] + δ
        if target < ε[1] || target > ε[end]
            continue
        elseif target == ε[end]
            W[e, end] = 1.0
            continue
        end
        j = searchsortedlast(ε, target)
        j = clamp(j, 1, N - 1)
        θ = (target - ε[j]) / (ε[j+1] - ε[j])
        W[e, j] = 1 - θ
        W[e, j+1] = θ
    end
    return W
end

"""
    build_conservative_shift_pair(ε, weights, δ)

Finite-volume translations on energy control volumes. The cells are centred
on the supplied nodes with interfaces at adjacent midpoints; endpoint cells
end at the first/last node. `weights` must be the corresponding trapezoidal
weights. A shifted cell receives its exact overlap with source cells, so this
is a positive, piecewise-constant conservative remap, **not** nodal linear
interpolation. For arbitrary displacement and nonuniform grids it satisfies
`Diagonal(weights)*plus == minus'*Diagonal(weights)` to roundoff. Values
outside the window are zero; `loss_plus/minus` expose fractional edge loss.
No periodic wrap, energy rounding, or row renormalisation is performed.
"""
function build_conservative_shift_pair(
    ε::AbstractVector{<:Real},
    weights::AbstractVector{<:Real},
    δ::Real,
)
    N = length(ε)
    N >= 2 || throw(ArgumentError("at least two energy nodes are required"))
    length(weights) == N || throw(DimensionMismatch("energy weights differ"))
    all(isfinite, ε) && all(>(0), diff(ε)) && isfinite(δ) ||
        throw(ArgumentError("invalid energy grid or displacement"))
    edges = vcat(Float64(ε[1]), (ε[1:(end-1)] .+ ε[2:end]) ./ 2, Float64(ε[end]))
    widths = diff(edges)
    all(isfinite, weights) && all(>(0), weights) ||
        throw(ArgumentError("energy weights must be finite and positive"))
    # Uniform quadrature weights and midpoint differences differ by a few
    # ulps of the *node* scale on long grids. Comparing only relative to one
    # tiny cell width incorrectly rejects legitimate fine energy meshes.
    coordinate_roundoff = 32eps(Float64) * max(maximum(abs, ε), ε[end]-ε[1])
    all(abs.(weights .- widths) .<= coordinate_roundoff .+ 32eps(Float64) .* widths) ||
        throw(ArgumentError("weights must equal trapezoidal control-volume widths"))
    plus = zeros(Float64, N, N)
    for i = 1:N, j = 1:N
        overlap = max(0.0, min(edges[i+1] + δ, edges[j+1]) - max(edges[i] + δ, edges[j]))
        plus[i, j] = overlap / weights[i]
    end
    minus = Matrix(Diagonal(1 ./ weights) * plus' * Diagonal(weights))
    return (
        plus = plus,
        minus = minus,
        loss_plus = max.(0.0, 1 .- vec(sum(plus; dims = 2))),
        loss_minus = max.(0.0, 1 .- vec(sum(minus; dims = 2))),
        weighted_adjoint_residual = norm(
            Diagonal(weights) * plus - minus' * Diagonal(weights),
        ),
        discretization = :finite_volume_piecewise_constant,
    )
end

"""
Apply a dense energy-shift matrix to an `(N_E,N_k,N_b,N_b)` field.

See [Field-periodic closure](@ref theory-periodicity) and
[the energy-shift equation](@ref eq-energy-shift).
"""
function apply_energy_shift(W::AbstractMatrix{<:Real}, X::AbstractArray{<:Number,4})
    NE, Nk, Nb, Nb2 = size(X)
    Nb == Nb2 || throw(DimensionMismatch("last two axes must be square"))
    size(W) == (NE, NE) || throw(DimensionMismatch("shift matrix has wrong shape"))
    Y = zeros(promote_type(eltype(W), eltype(X)), size(X))
    for m = 1:Nk, a = 1:Nb, b = 1:Nb
        Y[:, m, a, b] .= W * view(X, :, m, a, b)
    end
    return Y
end

function _empty_selfenergy(NE::Int, Nk::Int, Nb::Int)
    shape = (NE, Nk, Nb, Nb)
    return SelfEnergyFamily(
        zeros(ComplexF64, shape),
        zeros(ComplexF64, shape),
        zeros(ComplexF64, shape),
    )
end

function _embedding_component(T::AbstractMatrix, shifted::AbstractArray{<:Number,4})
    NE, Nk, Nb, _ = size(shifted)
    Σ = zeros(ComplexF64, NE, Nk, Nb, Nb)
    for e = 1:NE, m = 1:Nk
        Σ[e, m, :, :] .= T * Matrix(view(shifted, e, m, :, :)) * T'
    end
    return Σ
end

"""
    embedding_self_energy(problem, green)

Build the **legacy bulk-feedback approximation** for all
three Keldysh components.  Returns `(total, plus, minus)`; every family has
shape `(N_E,N_k,N_b,N_b)`.

The neighbouring full bulk Green function is not a cavity/surface Green
function. This closure is retained to reproduce old calculations; it must
not be identified as an exact nearest-neighbour embedding. See
`finite_chain_green` for an independent open-chain reference.

Implements [EQ-EMBED-001](@ref eq-embedding).
See [Field-periodic closure](@ref theory-periodicity) for the neighbouring
periods represented by the returned `plus` and `minus` families.
"""
function embedding_self_energy(problem::NEGFProblem, green::GreenState)
    bp = problem.basis
    arrays = ((green.Gᴿ, :Σᴿ), (green.Gˡ, :Σˡ), (green.Gᵍ, :Σᵍ))
    plus_values = Dict{Symbol,Array{ComplexF64,4}}()
    minus_values = Dict{Symbol,Array{ComplexF64,4}}()
    for (G, name) in arrays
        plus_values[name] = _embedding_component(bp.T₊, apply_energy_shift(problem.W₊ᴱᵖ, G))
        minus_values[name] =
            _embedding_component(bp.T₋, apply_energy_shift(problem.W₋ᴱᵖ, G))
    end
    plus = SelfEnergyFamily(plus_values[:Σᴿ], plus_values[:Σˡ], plus_values[:Σᵍ])
    minus = SelfEnergyFamily(minus_values[:Σᴿ], minus_values[:Σˡ], minus_values[:Σᵍ])
    total = SelfEnergyFamily(plus.Σᴿ + minus.Σᴿ, plus.Σˡ + minus.Σˡ, plus.Σᵍ + minus.Σᵍ)
    return total, plus, minus
end

"""
    finite_chain_green(energy, onsite, couplings; sigma_retarded,
                       sigma_lesser, sigma_greater=nothing, centre, return_full=false)

Independent finite open block-chain reference. `couplings[j]` is the actual
Hamiltonian block `H[j,j+1]`; field offsets belong in `onsite[j]` and are never
rounded to energy-grid nodes. Local scattering/reservoir self-energies are
block diagonal. Left/right *cavities* are eliminated separately by exact
Schur complements, including their lesser and greater injection terms.

No bulk Green function is reused as a surface. No hidden regulator is added:
a physical or mathematical broadening must be supplied through all three
self-energy components. If greater is omitted it is derived by the Keldysh
identity. Returns centre `retarded/lesser/greater`, the two embeddings, and
optionally independently inverted full matrices for small reference tests.
This finite domain requires a size-convergence study before describing an
infinite system, and assumes no inter-cell scattering self-energy.
"""
function finite_chain_green(
    energy::Real,
    onsite::AbstractVector{<:AbstractMatrix},
    couplings::AbstractVector{<:AbstractMatrix};
    sigma_retarded = [zeros(ComplexF64, size(h)) for h in onsite],
    sigma_lesser = [zeros(ComplexF64, size(h)) for h in onsite],
    sigma_greater = nothing,
    centre::Integer = cld(length(onsite), 2),
    return_full::Bool = false,
)
    count = length(onsite)
    count > 0 && 1 <= centre <= count || throw(ArgumentError("invalid chain centre"))
    length(couplings) == count - 1 || throw(DimensionMismatch("chain couplings differ"))
    Nb = size(first(onsite), 1)
    Nb > 0 || throw(ArgumentError("empty chain block"))
    all(
        h ->
            size(h) == (Nb, Nb) &&
            all(isfinite, h) &&
            norm(h-h') <= 128eps(Float64)*max(norm(h), 1.0),
        onsite,
    ) || throw(ArgumentError("onsite blocks must be equally sized and Hermitian"))
    local_h = [Matrix{ComplexF64}((h+h')/2) for h in onsite]
    all(t -> size(t) == (Nb, Nb), couplings) ||
        throw(DimensionMismatch("coupling block differs"))
    for sigma in (sigma_retarded, sigma_lesser)
        length(sigma) == count && all(s -> size(s) == (Nb, Nb), sigma) ||
            throw(DimensionMismatch("local self-energy blocks differ"))
    end
    greater =
        sigma_greater === nothing ?
        [sigma_retarded[j] - sigma_retarded[j]' + sigma_lesser[j] for j = 1:count] :
        sigma_greater
    length(greater) == count && all(s -> size(s) == (Nb, Nb), greater) ||
        throw(DimensionMismatch("greater self-energy blocks differ"))
    eye = Matrix{ComplexF64}(I, Nb, Nb)
    zero_block = zeros(ComplexF64, Nb, Nb)
    function cavity(indices, direction)
        R = copy(zero_block)
        L = copy(zero_block)
        G = copy(zero_block)
        for j in indices
            coupling =
                direction == :left ? (j == 1 ? zero_block : couplings[j-1]') :
                (j == count ? zero_block : couplings[j])
            Rnew = inv(
                energy * eye - local_h[j] - sigma_retarded[j] - coupling * R * coupling',
            )
            Lnew = Rnew * (sigma_lesser[j] + coupling * L * coupling') * Rnew'
            Gnew = Rnew * (greater[j] + coupling * G * coupling') * Rnew'
            R, L, G = Rnew, Lnew, Gnew
        end
        return (retarded = R, lesser = L, greater = G)
    end
    left = cavity(1:(centre-1), :left)
    right = cavity(count:-1:(centre+1), :right)
    Tleft = centre == 1 ? zero_block : couplings[centre-1]'
    Tright = centre == count ? zero_block : couplings[centre]
    embedding(T, g) = (
        retarded = T * g.retarded * T',
        lesser = T * g.lesser * T',
        greater = T * g.greater * T',
    )
    minus = embedding(Tleft, left)
    plus = embedding(Tright, right)
    R = inv(
        energy * eye - local_h[centre] - sigma_retarded[centre] - minus.retarded -
        plus.retarded,
    )
    L = R * (sigma_lesser[centre] + minus.lesser + plus.lesser) * R'
    G = R * (greater[centre] + minus.greater + plus.greater) * R'
    full = nothing
    if return_full
        D = zeros(ComplexF64, count * Nb, count * Nb)
        SL, SG = similar(D), similar(D)
        fill!(SL, 0)
        fill!(SG, 0)
        for j = 1:count
            rows = ((j-1)*Nb+1):(j*Nb)
            D[rows, rows] .= energy * eye - local_h[j] - sigma_retarded[j]
            SL[rows, rows] .= sigma_lesser[j]
            SG[rows, rows] .= greater[j]
            if j < count
                next = (j*Nb+1):((j+1)*Nb)
                D[rows, next] .= -couplings[j]
                D[next, rows] .= -couplings[j]'
            end
        end
        fullR = inv(D)
        full =
            (retarded = fullR, lesser = fullR * SL * fullR', greater = fullR * SG * fullR')
    end
    return (retarded = R, lesser = L, greater = G, plus = plus, minus = minus, full = full)
end

"""Retarded scalar nearest-neighbour surface Green function in the upper half-plane."""
function chain_surface_green(z::Complex, hopping::Number)
    imag(z) > 0 ||
        throw(ArgumentError("retarded argument must have positive imaginary part"))
    t = abs(hopping)
    t == 0 && return inv(z)
    root = sqrt(z - 2t) * sqrt(z + 2t)
    return 2 / (z + root)
end

"""
    cavity_embedding_self_energy(problem, local_scattering, h; periods=8)

Controlled open-chain cavity embedding for the current field-periodic local
scattering approximation. The central cell is surrounded by `periods` cells
on each side. Cell j has `h_j=h_0-j*Eperiod` and all three local self-energies
`Σ_j(E)=Σ_0(E+j*Eperiod)`, using the same noncyclic linear interpolation.
No previous bulk Green function is fed back as a surface and no extra
broadening is introduced. Returns `(total,plus,minus)` as SelfEnergyFamily.

Outside the sampled energy window local Σ is zero. Convergence in BOTH the
number of cells and the energy-window width is mandatory; finite-chain
reflection and neglected nonlocal scattering remain explicit model errors.
The function is a correctness/reference mode, not an optimized large-window
production implementation.
"""
function cavity_embedding_self_energy(
    problem::NEGFProblem,
    local_scattering::SelfEnergyFamily,
    h::AbstractArray{<:Number,3};
    periods::Int = 8,
    energy_shift::Symbol = :dense,
    worker_count::Int = 1,
)
    periods >= 1 ||
        throw(ArgumentError("cavity must include at least one neighbouring cell"))
    worker_count >= 1 || throw(ArgumentError("cavity worker_count must be positive"))
    energy_shift in (:dense, :sparse_plan, :conservative_pair) ||
        throw(ArgumentError("unsupported cavity energy shift"))
    NE, Nk, Nb, Nb2 = size(local_scattering.Σᴿ)
    Nb == Nb2 && size(h) == (Nk, Nb, Nb) ||
        throw(DimensionMismatch("cavity Hamiltonian differs"))
    size(local_scattering.Σˡ) == size(local_scattering.Σᵍ) == size(local_scattering.Σᴿ) ||
        throw(DimensionMismatch("cavity self-energy components differ"))
    ε = problem.grids.ε
    length(ε) == NE || throw(DimensionMismatch("cavity energy grid differs"))
    drop = _scaled_physics(problem.physical, problem.scales).Eᵖ
    plus, minus = _empty_selfenergy(NE, Nk, Nb), _empty_selfenergy(NE, Nk, Nb)
    eye = Matrix{ComplexF64}(I, Nb, Nb)
    zero_block = zeros(ComplexF64, Nb, Nb)
    edges = vcat(first(ε), (ε[1:(end-1)] .+ ε[2:end]) ./ 2, last(ε))
    function local_block(field, index, displacement, m)
        if energy_shift === :conservative_pair
            left, right = edges[index]+displacement, edges[index+1]+displacement
            width = problem.grids.wᴱ[index]
            value = copy(zero_block)
            first_source = clamp(searchsortedlast(edges, left), 1, NE)
            last_source = clamp(searchsortedfirst(edges, right), 1, NE)
            for j = first_source:last_source
                overlap = max(0.0, min(right, edges[j+1])-max(left, edges[j]))
                overlap == 0 && continue
                value .+= (overlap/width) .* view(field, j, m, :, :)
            end
            return value
        end
        target = ε[index]+displacement
        (target < ε[1] || target > ε[end]) && return copy(zero_block)
        target == ε[end] && return Matrix(view(field, NE, m, :, :))
        j = clamp(searchsortedlast(ε, target), 1, NE-1)
        fraction = (target-ε[j])/(ε[j+1]-ε[j])
        return (1-fraction) .* Matrix(view(field, j, m, :, :)) .+
               fraction .* Matrix(view(field, j+1, m, :, :))
    end
    couplings = [problem.basis.T₊ for _ = 1:(2periods)]
    function evaluate_block!(linear)
        e, m = (linear-1) % NE + 1, (linear-1) ÷ NE + 1
        onsite = [Matrix(view(h, m, :, :)) - j*drop*eye for j = (-periods):periods]
        sr = [local_block(local_scattering.Σᴿ, e, j*drop, m) for j = (-periods):periods]
        sl = [local_block(local_scattering.Σˡ, e, j*drop, m) for j = (-periods):periods]
        sg = [local_block(local_scattering.Σᵍ, e, j*drop, m) for j = (-periods):periods]
        result = finite_chain_green(
            ε[e],
            onsite,
            couplings;
            sigma_retarded = sr,
            sigma_lesser = sl,
            sigma_greater = sg,
            centre = periods+1,
        )
        plus.Σᴿ[e, m, :, :] .= result.plus.retarded
        plus.Σˡ[e, m, :, :] .= result.plus.lesser
        plus.Σᵍ[e, m, :, :] .= result.plus.greater
        minus.Σᴿ[e, m, :, :] .= result.minus.retarded
        minus.Σˡ[e, m, :, :] .= result.minus.lesser
        minus.Σᵍ[e, m, :, :] .= result.minus.greater
    end
    workers = min(worker_count, Base.Threads.nthreads(:default), NE*Nk)
    if workers == 1
        for linear = 1:(NE*Nk)
            evaluate_block!(linear)
        end
    else
        @sync for worker = 1:workers
            Base.Threads.@spawn for linear = worker:workers:(NE*Nk)
                evaluate_block!(linear)
            end
        end
    end
    total = SelfEnergyFamily(plus.Σᴿ + minus.Σᴿ, plus.Σˡ + minus.Σˡ, plus.Σᵍ + minus.Σᵍ)
    return total, plus, minus
end
