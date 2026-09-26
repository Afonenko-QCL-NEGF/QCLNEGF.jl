module Suite_T022
include("../support/common.jl")

@testset "Tutorial localized basis" begin
    p = reference_parameters()
    n = tutorial_numerics()
    s = ScaleSystem()
    g = build_grids(p, n, s)
    profiles = build_profiles(p, n, s, g)
    basis = build_localized_basis(p, n, s, g, profiles)
    @test size(basis.Φ) == (n.N_z, n.N_b)
    @test basis.Φ' * basis.Φ ≈ I atol=1e-11
    @test basis.T₋ ≈ basis.T₊' atol=1e-13
    @test norm(basis.H₀ - basis.H₀') / norm(basis.H₀) < 1e-12
    @test norm(basis.Z - basis.Z') / norm(basis.Z) < 1e-12
    @test norm(basis.M⁻¹ - basis.M⁻¹') / norm(basis.M⁻¹) < 1e-12

    unlocalized = build_localized_basis(p, n, s, g, profiles; localization = :none)
    @test basis.localization == :pzp
    @test unlocalized.localization == :none
    @test size(unlocalized.Φ) == (n.N_z, n.N_b)
    @test unlocalized.Φ' * unlocalized.Φ ≈ I atol=1e-11
    @test unlocalized.T₋ ≈ unlocalized.T₊' atol=1e-13
    @test norm(unlocalized.H₀ - unlocalized.H₀') / norm(unlocalized.H₀) < 1e-12
    @test_throws ArgumentError build_localized_basis(
        p,
        n,
        s,
        g,
        profiles;
        localization = :unknown,
    )
end

end # independent suite
