"""
    build_poisson_matrix(profiles, grids)

Form the full dense cyclically tridiagonal matrix of
`-d/dx(εᵣ du/dx)`.  The matrix is real symmetric, positive semidefinite,
and has rank `N_z-1` before the gauge row/column is added.

Implements [EQ-POISSON-001](@ref eq-poisson-matrix).
See [Periodic Poisson and the outer loop](@ref theory-outer-loop).
"""
function build_poisson_matrix(profiles::MaterialProfiles, grids::ModelGrids)
    N = length(grids.x)
    Δx = grids.wˣ[1]
    L = zeros(Float64, N, N)
    face = zeros(Float64, N)
    for g = 1:N
        gp = g == N ? 1 : g + 1
        face[g] = 2profiles.εᵣ[g] * profiles.εᵣ[gp] / (profiles.εᵣ[g] + profiles.εᵣ[gp])
    end
    for g = 1:N
        gm = g == 1 ? N : g - 1
        gp = g == N ? 1 : g + 1
        L[g, g] = (face[gm] + face[g]) / Δx^2
        L[g, gm] = -face[gm] / Δx^2
        L[g, gp] = -face[g] / Δx^2
    end
    return L
end

function _solve_periodic_poisson_with_diagnostics(
    profiles::MaterialProfiles,
    grids::ModelGrids,
    scales::ScaleSystem,
    n̄::AbstractVector{<:Real},
)
    N = length(grids.x)
    length(n̄) == N || throw(DimensionMismatch("density must have N_z entries"))
    Nᴰ²ᴰ = dot(grids.wˣ, profiles.Nᴰ)
    signed_charge = dot(grids.wˣ, profiles.Nᴰ .- n̄)
    absolute_charge = abs(signed_charge)
    neutrality = absolute_charge / Nᴰ²ᴰ
    L = build_poisson_matrix(profiles, grids)
    b = -scales.λ_P .* (profiles.Nᴰ .- n̄)
    sP = maximum(sum(abs, L; dims = 2))
    v = fill(inv(sqrt(N)), N)
    bordered = zeros(Float64, N + 1, N + 1)
    bordered[1:N, 1:N] .= L ./ sP
    bordered[1:N, end] .= v
    bordered[end, 1:N] .= v
    rhs = vcat(b ./ sP, 0.0)
    answer = bordered \ rhs
    Uᴴ = answer[1:N]
    ζ = answer[end]
    p_num = norm(L * Uᴴ - b)
    rhs_F = norm(b)
    p_floor = scales.λ_P * 1e-14
    p_den = rhs_F + p_floor
    rP = p_num / p_den
    return (;
        solution=(Uᴴ, ζ, rP, neutrality),
        diagnostics=(; residual_ratio=rP, residual_numerator_F=p_num, rhs_F,
                     additive_floor=p_floor, denominator=p_den, border_scale=sP,
                     zeta_signed=ζ, neutrality_signed=signed_charge,
                     neutrality_abs=absolute_charge, target_sheet=Nᴰ²ᴰ,
                     neutrality_ratio=neutrality),
    )
end

_solve_periodic_poisson_with_diagnostics(problem::NEGFProblem, n̄::AbstractVector{<:Real}) =
    _solve_periodic_poisson_with_diagnostics(problem.profiles, problem.grids, problem.scales, n̄)

"""
    solve_periodic_poisson(problem, n̄)

Solve the dimensionless bordered periodic Poisson system for the electronic
Hartree energy `Uᴴ/E₀=-φᴴ/φ₀`.  `n̄=L₀³n` must have `N_z` entries.  Returns
`(Uᴴ, ζ, r_P, r_neutral)`.

Implements [EQ-POISSON-002](@ref eq-bordered-poisson).
See [Periodic Poisson and the outer loop](@ref theory-outer-loop) for the
charge sign, periodic gauge, and neutrality condition.
"""
solve_periodic_poisson(
    profiles::MaterialProfiles, grids::ModelGrids, scales::ScaleSystem,
    n̄::AbstractVector{<:Real},
) = _solve_periodic_poisson_with_diagnostics(profiles, grids, scales, n̄).solution

solve_periodic_poisson(problem::NEGFProblem, n̄::AbstractVector{<:Real}) =
    _solve_periodic_poisson_with_diagnostics(problem, n̄).solution
