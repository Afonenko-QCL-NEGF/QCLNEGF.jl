function _face_inverse_mass(m_left::Real, m_right::Real)
    return 2.0 / (m_left + m_right)
end

"""
    build_bdd_hamiltonian(Eᶜ, mᶻᵣ, Δx, C̄; U=zeros, boundary=:dirichlet)

Return the full dense, real symmetric BenDaniel–Duke matrix.  Only its first
off-diagonals are nonzero.  This implementation deliberately forms the full
matrix before calling standard dense Hermitian eigensolvers.

Implements [EQ-BDD-001](@ref eq-bdd-matrix).
Physical background and matrix structure are given in
[BenDaniel–Duke and the fixed basis](@ref theory-basis).
"""
function build_bdd_hamiltonian(
    Eᶜ::AbstractVector{<:Real},
    mᶻᵣ::AbstractVector{<:Real},
    Δx::Real,
    C̄::Real;
    U::AbstractVector{<:Real} = zeros(length(Eᶜ)),
    boundary::Symbol = :dirichlet,
)
    N = length(Eᶜ)
    length(mᶻᵣ) == N == length(U) || throw(DimensionMismatch("profile lengths differ"))
    boundary in (:dirichlet, :periodic) || throw(ArgumentError("unknown boundary"))
    H = zeros(Float64, N, N)
    t = zeros(Float64, N - 1)
    for g = 1:(N-1)
        t[g] = C̄ * _face_inverse_mass(mᶻᵣ[g], mᶻᵣ[g+1]) / Δx^2
        H[g, g+1] = H[g+1, g] = -t[g]
    end
    t_left = C̄ / (mᶻᵣ[1] * Δx^2)
    t_right = C̄ / (mᶻᵣ[end] * Δx^2)
    for g = 1:N
        left = g == 1 ? t_left : t[g-1]
        right = g == N ? t_right : t[g]
        H[g, g] = Eᶜ[g] + U[g] + left + right
    end
    if boundary === :periodic
        tseam = C̄ * _face_inverse_mass(mᶻᵣ[end], mᶻᵣ[1]) / Δx^2
        H[1, 1] += tseam - t_left
        H[end, end] += tseam - t_right
        H[1, end] = H[end, 1] = -tseam
    end
    return H
end

function build_bdd_hamiltonian(
    profiles::MaterialProfiles,
    grids::ModelGrids,
    scales::ScaleSystem;
    U = zeros(length(grids.x)),
    boundary = :dirichlet,
)
    Δx = grids.wˣ[1]
    return build_bdd_hamiltonian(
        profiles.Eᶜ,
        profiles.mᶻᵣ,
        Δx,
        _C̄(scales);
        U = U,
        boundary = boundary,
    )
end

function _boundary_hopping(
    profiles::MaterialProfiles,
    grids::ModelGrids,
    scales::ScaleSystem,
)
    return _C̄(scales) * _face_inverse_mass(profiles.mᶻᵣ[end], profiles.mᶻᵣ[1]) /
           grids.wˣ[1]^2
end
