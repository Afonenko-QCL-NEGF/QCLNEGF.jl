module Suite_T106
include("../support/common.jl")

@testset "reference design units, geometry and quadratures" begin
    p = reference_parameters()
    n = baseline_numerics()
    s = ScaleSystem()
    g = build_grids(p, n, s)
    profiles = build_profiles(p, n, s, g)

    @test ustrip(u"nm", uconvert(u"nm", sum(layer.d for layer in p.layers))) ≈ 29.61
    @test length(g.x) == 198
    @test length(g.ε) == 2401
    @test length(g.trusted_energy) == 2401
    @test any(g.trusted_energy)
    @test sum(g.wᵏ) ≈ (6.0^2) / (4π)
    @test sum(g.wᴱ) ≈ 6.0
    @test sum(g.wᑫᶻ) ≈ 200.0
    @test count(!iszero, profiles.Nᴰ) == 19
    @test findall(!iszero, profiles.Nᴰ) == collect(145:163)
    target = QCLNEGF._scaled_physics(p, s).Nᴰ²ᴰ
    @test dot(g.wˣ, profiles.Nᴰ) ≈ target rtol=2e-15
    physical_ND = maximum(profiles.Nᴰ) / s.L₀_m^3
    @test physical_ND ≈ 1.5837466005e23 rtol=2e-9

    @test scaled_value(100u"meV", s, :energy) ≈ 1.0
    @test QCLNEGF.C₀ ==
          uconvert(u"eV*nm^2", QCLNEGF.CODATA.ħ^2 / (2 * QCLNEGF.CODATA.m₀))
    @test QCLNEGF.CODATA.ħₑᵥ == uconvert(u"eV*s", QCLNEGF.CODATA.ħ)
    @test QCLNEGF.CODATA.kᴮₑᵥ == uconvert(u"eV/K", QCLNEGF.CODATA.kᴮ)
    @test ustrip(u"meV", uconvert(u"meV", physical_value(1.0, s, :energy))) ≈ 100.0
    @test_throws Unitful.DimensionError reference_parameters(F_bias = 1u"K")
    Ep = p.F_bias * sum(layer.d for layer in p.layers)
    @test ustrip(u"mV", uconvert(u"mV", Ep)) ≈ 56.0 rtol=2e-15
    @test_throws ArgumentError ScaleSystem(E₀ = -0.1u"eV")
    @test_throws ArgumentError reference_parameters(Tᴸ = 0u"K")
    @test_throws ArgumentError reference_parameters(ΔV_alloy = 0.6u"eV")
    palloy = reference_parameters(ΔV_alloy = 0.6u"eV", Ω₀ = 4.5e-29u"m^3")
    @test palloy.ΔV_alloy == 0.6u"eV"
    @test palloy.Ω₀ == 4.5e-29u"m^3"
    @test_throws ArgumentError QCLNEGF._check_options(
        SolverOptions(tolerances = SolverTolerances(r_D = NaN)),
    )
    @test_throws ArgumentError QCLNEGF._check_options(
        SolverOptions(energy_tail_window_fraction = 1.0),
    )
end

end # independent suite
