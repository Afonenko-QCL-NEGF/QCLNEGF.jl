"""
    build_poisson_matrix(profiles, grids)

Построить полную плотную циклическую трёхдиагональную матрицу оператора
`-d/dx(εᵣ du/dx)`. До добавления окаймляющей строки и столбца
матрица вещественная, симметричная и положительно полуопределённая
с рангом `N_z-1`.

Реализует [EQ-POISSON-001](@ref eq-poisson-matrix).
См. [периодическое уравнение Пуассона и внешний цикл](@ref theory-outer-loop).
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

"""
    solve_periodic_poisson(problem, n̄)

Решить безразмерную окаймлённую периодическую систему Пуассона
для энергии Хартри электрона `Uᴴ/E₀=-φᴴ/φ₀`. Вектор `n̄=L₀³n`
содержит `N_z` элементов. Возвращает `(Uᴴ, ζ, r_P, r_neutral)`.

Реализует [EQ-POISSON-002](@ref eq-bordered-poisson). См. знак заряда,
фиксацию постоянной и нейтральность в [главе о внешнем цикле](@ref theory-outer-loop).
"""
function solve_periodic_poisson(
    profiles::MaterialProfiles,
    grids::ModelGrids,
    scales::ScaleSystem,
    n̄::AbstractVector{<:Real},
)
    N = length(grids.x)
    length(n̄) == N || throw(DimensionMismatch("density must have N_z entries"))
    Nᴰ²ᴰ = dot(grids.wˣ, profiles.Nᴰ)
    neutrality = abs(dot(grids.wˣ, profiles.Nᴰ .- n̄)) / Nᴰ²ᴰ
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
    rP = norm(L * Uᴴ - b) / (norm(b) + scales.λ_P * 1e-14)
    return Uᴴ, ζ, rP, neutrality
end

solve_periodic_poisson(problem::NEGFProblem, n̄::AbstractVector{<:Real}) =
    solve_periodic_poisson(problem.profiles, problem.grids, problem.scales, n̄)
