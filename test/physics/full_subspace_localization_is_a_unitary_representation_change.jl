module Suite_T021
include("../support/common.jl")

@testset "Full-subspace localization is a unitary representation change" begin
    p = reference_parameters()
    template = tutorial_numerics()
    n = NumericalParameters(
        N_z = template.N_z,
        N_b = template.N_z,
        P_basis = template.P_basis,
        E_min = template.E_min,
        E_max = template.E_max,
        N_E = template.N_E,
        M_E = template.M_E,
        k_max = template.k_max,
        N_k = template.N_k,
        N_φ = template.N_φ,
        qz_max = template.qz_max,
        N_qz = template.N_qz,
        η_seed = template.η_seed,
    )
    s = ScaleSystem()
    g = build_grids(p, n, s)
    profiles = build_profiles(p, n, s, g)
    pzp = build_basis(p, n, s, g, profiles; localization = :pzp)
    energy = build_basis(p, n, s, g, profiles; localization = :none)
    U = pzp.Φ' * energy.Φ

    @test U' * U ≈ I atol=5e-11
    @test energy.H₀ ≈ U' * pzp.H₀ * U rtol=5e-11 atol=5e-11
    @test energy.Z ≈ U' * pzp.Z * U rtol=5e-11 atol=5e-11
    @test energy.M⁻¹ ≈ U' * pzp.M⁻¹ * U rtol=5e-11 atol=5e-11
    @test energy.T₊ ≈ U' * pzp.T₊ * U rtol=5e-11 atol=5e-11
end

end # independent suite
