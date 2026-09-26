module Suite_T073
include("../support/common.jl")

@testset "Manufactured periodic Poisson solution" begin
    p = reference_parameters()
    n = tutorial_numerics()
    s = ScaleSystem()
    g = build_grids(p, n, s)
    profiles = build_profiles(p, n, s, g)
    L = build_poisson_matrix(profiles, g)
    Utest = 0.1 .* cos.(2π .* g.x ./ QCLNEGF._scaled_physics(p, s).Lp)
    Utest .-= sum(g.wˣ .* Utest) / sum(g.wˣ)
    nbar = profiles.Nᴰ .+ (L * Utest) ./ s.λ_P
    U, ζ, rP, rneutral = solve_periodic_poisson(profiles, g, s, nbar)
    @test U ≈ Utest atol=2e-12
    @test abs(ζ) < 1e-12
    @test rP < 1e-11
    @test rneutral < 1e-12
end

end # independent suite
