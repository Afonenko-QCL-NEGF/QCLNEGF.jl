module ElectronElectronTests
using Test
using LinearAlgebra
using Unitful
using QCLNEGF

const Models = QCLNEGF.PhysicalModelExtensions
const Operators = QCLNEGF.QCLReferenceOperators

@testset "SPPA pole: primary analytical limits and parameter contract" begin
    disabled = ElectronElectronOptions()
    @test disabled.mode === :none
    @test all(
        isnothing,
        (
            disabled.transfer_wavenumber,
            disabled.electron_temperature,
            disabled.effective_mass_ratio,
            disabled.carrier_density_si,
        ),
    )
    @test_throws ArgumentError ElectronElectronOptions(
        :none,
        2,
        1u"m^-1",
        1u"K",
        1.0,
        1.0,
        true,
    )
    @test_throws ArgumentError ElectronElectronOptions(
        :sppa_single_q,
        2,
        nothing,
        nothing,
        nothing,
        nothing,
        true,
    )
    mass = 0.067 * CODATA.m₀
    screening = Models.DebyeScreening3D(1e20u"m^-3", 70u"K", 12.9)
    pole = Models.single_plasmon_pole(screening, 1e3u"m^-1", mass, 70u"K")
    plasma_oracle = CODATA.ħ * sqrt(CODATA.e^2 * 1e20u"m^-3" / (CODATA.ε₀ * 12.9 * mass))
    @test pole.plasma_energy ≈ plasma_oracle
    @test pole.pole_energy ≈ plasma_oracle rtol = 1e-6
    @test pole.correlation_weight ≈ plasma_oracle / 2 rtol = 1e-6
    s2 = Models.ThomasFermiScreening2D(4.5e14u"m^-2", 70u"K", mass, 2, 12.9)
    p1 = Models.single_plasmon_pole(s2, 1e5u"m^-1", mass, 70u"K")
    p4 = Models.single_plasmon_pole(s2, 4e5u"m^-1", mass, 70u"K")
    @test p4.plasma_energy ≈ 2p1.plasma_energy
    @test p1.pole_energy > p1.plasma_energy
    @test_throws ArgumentError Models.single_plasmon_pole(s2, 1e5u"m^-1", mass, 100u"K")
    @test_throws ArgumentError ElectronElectronOptions(mode = :sppa_single_q)
    @test_throws Unitful.DimensionError ElectronElectronOptions(
        mode = :sppa_single_q,
        screening_dimension = 2,
        transfer_wavenumber = 0.1u"nm^-1",
        electron_temperature = 70u"K",
        effective_mass_ratio = 0.067,
        carrier_density = 1e20u"m^-3",
    )
end

@testset "Actual SPPA SCBA family and independent subband quadrature" begin
    ee = ElectronElectronOptions(
        mode = :sppa_single_q,
        screening_dimension = 2,
        transfer_wavenumber = 0.1u"nm^-1",
        electron_temperature = 70u"K",
        effective_mass_ratio = 0.067,
        carrier_density = 4.5e14u"m^-2",
    )
    physical = reference_parameters(Tᴸ = 70u"K")
    scattering = ScatteringOptions(
        LO = false,
        acoustic = false,
        impurity = false,
        IFR = false,
        alloy = false,
    )
    problem = build_problem(
        physical = physical,
        numerical = tutorial_numerics(),
        scattering = scattering,
        physical_models = PhysicalModelOptions(electron_electron = ee),
        validate_static = false,
    )
    @test :electron_electron in problem.kernels.enabled
    @test problem.kernels.qᴷ[:electron_electron] > 0
    qbar = ustrip(u"m^-1", ee.transfer_wavenumber) * problem.scales.L₀_m
    form = Operators.sppa_form_factor(problem.grids, problem.basis, qbar)
    χ, x, w = problem.basis.χ, problem.grids.x, problem.grids.wˣ
    # Direct double quadrature checks conjugation and index ordering of the
    # matrix-based implementation on an off-diagonal subband element.
    a, c, d, b = 1, 2, 2, 1
    direct = sum(
        w[g] *
        w[h] *
        conj(χ[g, a]) *
        χ[g, c] *
        exp(-qbar*abs(x[g]-x[h])) *
        conj(χ[h, d]) *
        χ[h, b] for g in eachindex(x), h in eachindex(x)
    )
    @test form[a, c, d, b] ≈ direct rtol = 1e-12
    nb = problem.numerical.N_b
    covariance = [form[a, c, d, b] for a = 1:nb, c = 1:nb, b = 1:nb, d = 1:nb]
    covariance = reshape(covariance, nb^2, nb^2)
    @test covariance ≈ covariance'
    @test minimum(eigvals(Hermitian(covariance))) >= -1e-12 * opnorm(covariance)

    ne, nk = problem.numerical.N_E, problem.numerical.N_k
    hamiltonian = zeros(ComplexF64, nk, nb, nb)
    sigma = zeros(ComplexF64, ne, nk, nb, nb)
    for m = 1:nk, a = 1:nb
        hamiltonian[m, a, a] = 0.2a
        sigma[:, m, a, a] .= -0.1im
    end
    gr, condition, scaling =
        retarded_green(problem.grids.ε, hamiltonian, sigma; return_scale = true)
    spectral = spectral_function(gr)
    lesser = zeros(ComplexF64, size(gr))
    thermal_energy = ustrip(u"eV", CODATA.kᴮₑᵥ * 70u"K") / problem.scales.E₀_eV
    for e = 1:ne, m = 1:nk, a = 1:nb
        occupation = 1 / (1 + exp((problem.grids.ε[e]-0.5) / thermal_energy))
        lesser[e, m, a, a] = im * occupation * spectral[e, m, a, a]
    end
    greater = greater_green(gr, lesser)
    green = GreenState(gr, lesser, greater, spectral, condition, scaling)
    parts = Operators.sppa_components(problem, green)
    @test size(parts.lesser) == size(gr)
    @test norm(parts.lesser) > 0
    @test norm(parts.greater) > 0
    @test parts.exchange ≈ parts.exchange'
    @test maximum(eigvals(Hermitian(parts.exchange))) <= 1e-10
    for e = 1:ne, m = 1:nk
        @test minimum(eigvals(Hermitian(-im * parts.lesser[e, m, :, :]))) >= -1e-11
        @test minimum(eigvals(Hermitian(im * parts.greater[e, m, :, :]))) >= -1e-11
    end
    family = Operators.sppa_self_energy(problem, green)
    @test family isa SelfEnergyFamily
    @test all(isfinite, family.Σᴿ)
    # Causality identity remains independent of the static Hermitian Fock term.
    for e = 1:ne, m = 1:nk
        @test family.Σᴿ[e, m, :, :] - family.Σᴿ[e, m, :, :]' ≈
              family.Σᵍ[e, m, :, :] - family.Σˡ[e, m, :, :] atol = 1e-10
    end
    diagnostics = Operators.electron_electron_collision_diagnostics(problem, green, family)
    @test isfinite(diagnostics.energy_residual)
    @test !diagnostics.conserving_GW_certificate

    # These comparisons verify dispatch wiring for the shared SPPA model;
    # independent scientific anchors are the pole and quadrature tests above.
    reference_candidates = Operators._scattering_candidate(problem, green)
    @test reference_candidates[:electron_electron].Σˡ ≈ family.Σˡ
    options = ProductionOptions(algorithms = AlgorithmOptions(hilbert = :direct))
    cache = build_production_cache(problem; options)
    production_candidates = QCLNEGF.QCLNumerics._scattering_candidate_production(
        problem,
        green,
        cache,
        options,
    )
    @test production_candidates[:electron_electron].Σᴿ ≈ family.Σᴿ
    @test production_candidates[:electron_electron].Σˡ ≈ family.Σˡ

    conservative = NEGFProblem(
        problem.physical,
        problem.numerical,
        problem.scattering,
        problem.scales,
        problem.grids,
        problem.profiles,
        problem.basis,
        problem.kernels,
        problem.W₊ᴱᵖ,
        problem.W₋ᴱᵖ,
        problem.W₊ᴸᴼ,
        problem.W₋ᴸᴼ,
        :finite_volume_piecewise_constant,
        problem.models,
    )
    conserved_parts = Operators.sppa_components(conservative, green)
    collision_family = SelfEnergyFamily(
        zeros(ComplexF64, size(gr)),
        conserved_parts.lesser,
        conserved_parts.greater,
    )
    conserved = Operators.electron_electron_collision_diagnostics(
        conservative,
        green,
        collision_family,
    )
    # Weighted-adjoint translations conserve number even when energy transfer
    # to the declared equilibrium plasmon closure is nonzero.
    @test conserved.charge_residual < 1e-12

    # The disabled option keeps the exact existing kernel object and model.
    @test Operators.add_sppa_kernel(
        physical,
        problem.scales,
        problem.grids,
        problem.basis,
        problem.kernels,
        ElectronElectronOptions(),
    ) === problem.kernels
end
end
