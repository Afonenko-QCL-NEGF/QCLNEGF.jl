_matrix_block(X, e, m) = Matrix(view(X, e, m, :, :))

"""
    retarded_green(ε, h, Σᴿ; η=0, return_scale=false)

Directly invert one full dense Dyson matrix for every `(E_e,k_m)`.  A scalar
row-sum rescaling is applied before `inv`; it is algebraically exact and does
not change the condition number.  The returned `Gᴿ` is dimensionless and has
shape `(N_E,N_k,N_b,N_b)`.
With `return_scale=true`, the third result is the dimensionless row-sum
equilibration scale ``s_{em}``; otherwise the stable public return remains
`(Gᴿ, condition_number)`.

Implements [EQ-DYSON-001](@ref eq-dyson-retarded).
See [Green functions and density](@ref theory-greens) for the index and sign
conventions.
"""
function retarded_green(
    ε::AbstractVector{<:Real},
    h::AbstractArray{<:Number,3},
    Σᴿ::AbstractArray{<:Number,4};
    η::Real = 0.0,
    return_scale::Bool = false,
)
    NE, Nk, Nb, Nb2 = size(Σᴿ)
    Nb == Nb2 || throw(DimensionMismatch("Σᴿ blocks must be square"))
    size(h) == (Nk, Nb, Nb) || throw(DimensionMismatch("h has wrong shape"))
    length(ε) == NE || throw(DimensionMismatch("energy axis differs"))
    Gᴿ = zeros(ComplexF64, NE, Nk, Nb, Nb)
    κD = zeros(Float64, NE, Nk)
    scaleD = zeros(Float64, NE, Nk)
    I_b = Matrix{ComplexF64}(I, Nb, Nb)
    for e = 1:NE, m = 1:Nk
        D = (ε[e] + im * η) * I_b - Matrix(view(h, m, :, :)) - _matrix_block(Σᴿ, e, m)
        scale = maximum(sum(abs, D; dims = 2))
        scale > 0 || throw(SingularException(0))
        Gᴿ[e, m, :, :] .= inv(D / scale) / scale
        κD[e, m] = cond(D)
        scaleD[e, m] = scale
    end
    return return_scale ? (Gᴿ, κD, scaleD) : (Gᴿ, κD)
end

"""
Return `A=i(Gᴿ-Gᴬ)` for every energy and momentum block.

See [Green functions and density](@ref theory-greens) and
[the Keldysh identities](@ref eq-keldysh).
"""
function spectral_function(Gᴿ::AbstractArray{<:Number,4})
    A = similar(ComplexF64.(Gᴿ))
    for e in axes(Gᴿ, 1), m in axes(Gᴿ, 2)
        GR = _matrix_block(Gᴿ, e, m)
        A[e, m, :, :] .= im .* (GR - GR')
    end
    return A
end

"""
    keldysh_green(Gᴿ, Σˡ)

Evaluate the full dense Keldysh product `Gˡ=GᴿΣˡGᴬ` independently at every
`(E_e,k_m)`.  Implements [EQ-KELDYSH-001](@ref eq-keldysh).

See [Green functions and density](@ref theory-greens).
"""
function keldysh_green(Gᴿ::AbstractArray{<:Number,4}, Σˡ::AbstractArray{<:Number,4})
    size(Gᴿ) == size(Σˡ) || throw(DimensionMismatch("Gᴿ and Σˡ differ"))
    Gˡ = zeros(ComplexF64, size(Gᴿ))
    for e in axes(Gᴿ, 1), m in axes(Gᴿ, 2)
        GR = _matrix_block(Gᴿ, e, m)
        Gˡ[e, m, :, :] .= GR * _matrix_block(Σˡ, e, m) * GR'
    end
    return Gˡ
end

"""
Construct `Gᵍ=Gᴿ-Gᴬ+Gˡ` without storing an independent advanced field.

See [Green functions and density](@ref theory-greens) and
[the Keldysh identities](@ref eq-keldysh).
"""
function greater_green(Gᴿ::AbstractArray{<:Number,4}, Gˡ::AbstractArray{<:Number,4})
    size(Gᴿ) == size(Gˡ) || throw(DimensionMismatch("G arrays differ"))
    Gᵍ = zeros(ComplexF64, size(Gᴿ))
    for e in axes(Gᴿ, 1), m in axes(Gᴿ, 2)
        GR = _matrix_block(Gᴿ, e, m)
        Gᵍ[e, m, :, :] .= GR - GR' + _matrix_block(Gˡ, e, m)
    end
    return Gᵍ
end

"""
Dimensionless sheet-number functional of an `(E,k,a,b)` lesser field.
Implements [EQ-NUMBER-001](@ref eq-number-normalization).

See [Green functions and density](@ref theory-greens).
"""
function number_functional(Gˡ::AbstractArray{<:Number,4}, grids::ModelGrids, g_s::Integer)
    NE, Nk, _, _ = size(Gˡ)
    value = 0.0 + 0.0im
    for e = 1:NE, m = 1:Nk
        value += grids.wᴱ[e] * grids.wᵏ[m] * tr(_matrix_block(Gˡ, e, m))
    end
    value *= -im * g_s / (2π)
    abs(imag(value)) ≤ _NUMBER_FUNCTIONAL_REALNESS_LIMIT * max(abs(real(value)), 1.0) ||
        throw(DomainError(value, "number functional is not real"))
    return real(value)
end

"""
    normalize_lesser(Gˡ, target, grids, g_s)

Legacy unilateral number normalization for explicit comparison only. Return
`(Gˡ_N, λ)`, where `λ=N[Gˡ]/target`. This operation does not preserve the
upper Pauli bound and is no longer the default SCBA iteration. A physical
fixed point still requires `λ→1`; the default uses `_normalize_keldysh_pair!`.
Implements [EQ-NUMBER-001](@ref eq-number-normalization).
See [Green functions and density](@ref theory-greens) for why the scalar
normalization eigenvalue must tend to one at the physical fixed point.
"""
function normalize_lesser(
    Gˡ::AbstractArray{<:Number,4},
    target::Real,
    grids::ModelGrids,
    g_s::Integer,
)
    isfinite(target) && target > 0 ||
        throw(ArgumentError("target charge must be finite and positive"))
    λ = number_functional(Gˡ, grids, g_s) / target
    isfinite(λ) && λ > 0 || throw(DomainError(λ, "kinetic eigenvalue λ must be positive"))
    return ComplexF64.(Gˡ ./ λ), λ
end

function _green_state(Gᴿ, Gˡ, κD, scaleD)
    Gᵍ = greater_green(Gᴿ, Gˡ)
    A = spectral_function(Gᴿ)
    return GreenState(ComplexF64.(Gᴿ), ComplexF64.(Gˡ), Gᵍ, A, κD, scaleD)
end

"""Coefficients for a fixed-number, Pauli-preserving Keldysh update.

Let `N` and `H` be the integrals of the positive occupied and empty matrices
`Gn=-iG<` and `Gp=iG>`. For excess electrons, transfer `(1-target/N)Gn`
to `Gp`; for missing electrons, transfer `(target-N)/H` of `Gp` to `Gn`.
Both maps are convex and preserve `Gn+Gp`, so they preserve the Pauli cone
when the raw candidate belongs to it. No eigenvalues are clipped. The map
becomes the identity exactly when the raw number equals the target; the
separate raw `lambda=N/target` gate remains mandatory. This is an iteration
constraint, not a finite-bias Fermi distribution or a physical reservoir.
"""
function _paired_number_coefficients(number::Real, holes::Real, target::Real)
    isfinite(target) && target > 0 ||
        throw(ArgumentError("target charge must be finite and positive"))
    isfinite(number) && number > 0 ||
        throw(DomainError(number, "raw occupied number must be finite and positive"))
    isfinite(holes) && holes >= 0 ||
        throw(DomainError(holes, "raw empty-state capacity must be finite and nonnegative"))
    target <= number + holes || throw(
        DomainError(
            (;
                target,
                raw_occupied = number,
                raw_unoccupied = holes,
                capacity = number+holes,
            ),
            "target charge exceeds the represented Pauli capacity (solver units); refine the domain/basis",
        ),
    )
    if number >= target
        return (;
            occupied_fraction = target/number,
            empty_fraction = 0.0,
            lambda = number/target,
        )
    end
    return (;
        occupied_fraction = 1.0,
        empty_fraction = (target-number)/holes,
        lambda = number/target,
    )
end

"""Apply the paired number constraint in place to independently evaluated G</G>."""
function _normalize_keldysh_pair!(Gˡ, Gᵍ, number::Real, holes::Real, target::Real)
    size(Gˡ) == size(Gᵍ) || throw(DimensionMismatch("Keldysh arrays differ"))
    weights = _paired_number_coefficients(number, holes, target)
    a, c = weights.occupied_fraction, weights.empty_fraction
    for index in eachindex(Gˡ, Gᵍ)
        lesser, greater = Gˡ[index], Gᵍ[index]
        # G< = i Gn and G> = -i Gp: transfers therefore have minus signs.
        Gˡ[index] = a*lesser - c*greater
        Gᵍ[index] = (1-c)*greater - (1-a)*lesser
    end
    return weights.lambda
end
