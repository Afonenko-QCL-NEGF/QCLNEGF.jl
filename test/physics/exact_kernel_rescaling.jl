module Suite_T048
include("../support/common.jl")

@testset "Exact kernel rescaling" begin
    p = reference_parameters()
    n = tutorial_numerics()
    s = ScaleSystem()
    grids = build_grids(p, n, s)
    profiles = build_profiles(p, n, s, grids)
    basis = build_localized_basis(p, n, s, grids, profiles)
    raw = acoustic_kernel(p, s, grids, basis)
    set = build_kernels(
        p,
        n,
        ScatteringOptions(
            LO = false,
            acoustic = true,
            impurity = false,
            IFR = false,
            alloy = false,
        ),
        s,
        grids,
        profiles,
        basis,
    )
    @test set.enabled == [:acoustic]
    @test maximum(abs, set.K[:acoustic]) ≈ 1.0 rtol=4eps(Float64)
    @test set.qᴷ[:acoustic] .* set.K[:acoustic] ≈ raw rtol=4eps(Float64)
    rH, rPSD = QCLNEGF._kernel_residual(set.K[:acoustic])
    @test rH < 1e-12
    @test rPSD < 1e-12

    full = build_kernels(p, n, default_scattering(), s, grids, profiles, basis)
    @test Set(full.enabled) == Set((:LO, :acoustic, :impurity, :IFR))
    for name in full.enabled
        rH, rPSD = QCLNEGF._kernel_residual(full.K[name])
        @test rH < 1e-10
        @test rPSD < 1e-10
        @test full.qᴷ[name] > 0
    end
    hzero = argmin(abs.(grids.qᶻ))
    @test full.Fᴸᴼ[hzero, :, :] ≈ I atol=5e-12

    # The two LO Keldysh components use opposite emission/absorption shifts.
    # Use energy-dependent sentinels so an accidental sign swap is observable.
    sp = QCLNEGF._scaled_physics(p, s)
    WplusEp = build_shift_matrix(grids.ε, +sp.Eᵖ)
    WminusEp = build_shift_matrix(grids.ε, -sp.Eᵖ)
    WplusLO = build_shift_matrix(grids.ε, +sp.ħωᴸᴼ)
    WminusLO = build_shift_matrix(grids.ε, -sp.ħωᴸᴼ)
    problem = NEGFProblem(
        p,
        n,
        default_scattering(),
        s,
        grids,
        profiles,
        basis,
        full,
        WplusEp,
        WminusEp,
        WplusLO,
        WminusLO,
    )
    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    Gless = zeros(ComplexF64, shape)
    Ggreater = zeros(ComplexF64, shape)
    for e = 1:n.N_E, m = 1:n.N_k, a = 1:n.N_b, b = 1:n.N_b
        Gless[e, m, a, b] = im * (e + 0.1m + 0.01a + 0.001b)
        Ggreater[e, m, a, b] = -im * (2e + 0.2m + 0.02a + 0.002b)
    end
    green = GreenState(
        zeros(ComplexF64, shape),
        Gless,
        Ggreater,
        zeros(ComplexF64, shape),
        ones(n.N_E, n.N_k),
        ones(n.N_E, n.N_k),
    )
    KLO = full.K[:LO]
    qK = full.qᴷ[:LO]
    Σless, Σgreater = QCLNEGF._lo_contraction(KLO, green, problem, qK)
    NLO = inv(exp(ustrip(u"eV", p.ħωᴸᴼ) / ustrip(u"eV", CODATA.kᴮₑᵥ * p.Tᴸᴼ)) - 1)
    expected_less =
        (NLO + 1) .* QCLNEGF._static_contraction(
            KLO,
            apply_energy_shift(WplusLO, Gless),
            grids.wᵏ,
            qK,
        ) .+
        NLO .* QCLNEGF._static_contraction(
            KLO,
            apply_energy_shift(WminusLO, Gless),
            grids.wᵏ,
            qK,
        )
    expected_greater =
        (NLO + 1) .* QCLNEGF._static_contraction(
            KLO,
            apply_energy_shift(WminusLO, Ggreater),
            grids.wᵏ,
            qK,
        ) .+
        NLO .* QCLNEGF._static_contraction(
            KLO,
            apply_energy_shift(WplusLO, Ggreater),
            grids.wᵏ,
            qK,
        )
    @test Σless ≈ expected_less rtol=2e-14
    @test Σgreater ≈ expected_greater rtol=2e-14

    palloy = reference_parameters(ΔV_alloy = 0.6u"eV", Ω₀ = 4.5e-29u"m^3")
    Kalloy = alloy_kernel(palloy, s, grids, profiles, basis)
    @test maximum(abs, Kalloy) > 0
    rHalloy, rPSDalloy = QCLNEGF._kernel_residual(Kalloy)
    @test rHalloy < 1e-10
    @test rPSDalloy < 1e-10
end

end # independent suite
