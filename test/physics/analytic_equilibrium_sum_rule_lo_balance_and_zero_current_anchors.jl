module Suite_T008
include("../support/common.jl")

@testset "Analytic equilibrium, sum-rule, LO balance, and zero-current anchors" begin
    @testset "A-EQUIL-01 scalar KMS/FDT" begin
        energy = collect(range(-0.40, 0.40; length = 4001))
        level = 0.03
        broadening = 0.07
        chemical_potential = -0.02
        kBT = 0.025

        hamiltonian = reshape(ComplexF64[level], 1, 1, 1)
        sigma_retarded = fill(-0.5im * broadening, length(energy), 1, 1, 1)
        green_retarded, _ = retarded_green(energy, hamiltonian, sigma_retarded)
        spectral = spectral_function(green_retarded)
        occupation = @. inv(exp((energy - chemical_potential) / kBT) + 1)
        green_lesser = similar(green_retarded)
        green_lesser[:, 1, 1, 1] .= im .* occupation .* spectral[:, 1, 1, 1]
        green_greater = greater_green(green_retarded, green_lesser)

        expected_greater = @. -im * (1 - occupation) * spectral[:, 1, 1, 1]
        @test green_greater[:, 1, 1, 1] ≈ expected_greater rtol=3e-13
        interior = findall(@. (occupation > 1e-7) & (occupation < 1 - 1e-7))
        kms_ratio =
            imag.(green_greater[interior, 1, 1, 1]) ./
            imag.(green_lesser[interior, 1, 1, 1])
        @test kms_ratio ≈ -exp.((energy[interior] .- chemical_potential) ./ kBT) rtol=2e-12
    end

    @testset "A-SUMRULE-01 finite-window Lorentzian" begin
        energy = collect(range(-10.0, 10.0; length = 40001))
        level = 0.17
        broadening = 0.08
        hamiltonian = reshape(ComplexF64[level], 1, 1, 1)
        sigma_retarded = fill(-0.5im * broadening, length(energy), 1, 1, 1)
        green_retarded, _ = retarded_green(energy, hamiltonian, sigma_retarded)
        spectral = real.(spectral_function(green_retarded)[:, 1, 1, 1])
        spacing = energy[2] - energy[1]
        weights = fill(spacing, length(energy))
        weights[[1, end]] ./= 2
        numerical_window_weight = dot(weights, spectral) / (2pi)
        exact_window_weight =
            (
                atan(2 * (last(energy) - level) / broadening) -
                atan(2 * (first(energy) - level) / broadening)
            ) / pi
        @test numerical_window_weight ≈ exact_window_weight atol=1e-5
        @test 0 < numerical_window_weight < 1
    end

    @testset "A-LO-DB-01 Bose/Fermi detailed balance" begin
        temperature = 200.0
        phonon_energy_eV = 0.0367
        kBT_eV = Float64(ustrip(u"eV", CODATA.kᴮₑᵥ * temperature * u"K"))
        bose = inv(exp(phonon_energy_eV / kBT_eV) - 1)
        @test (bose + 1) / bose ≈ exp(phonon_energy_eV / kBT_eV) rtol=2e-15

        chemical_potential = 0.01
        fermi(value) = inv(exp((value - chemical_potential) / kBT_eV) + 1)
        for lower_energy in (-0.05, 0.0, 0.05)
            upper_energy = lower_energy + phonon_energy_eV
            downward = (bose + 1) * fermi(upper_energy) * (1 - fermi(lower_energy))
            upward = bose * fermi(lower_energy) * (1 - fermi(upper_energy))
            @test downward ≈ upward rtol=3e-14
        end
    end

    @testset "A-CURRENT-01 zero-bias equilibrium collision current" begin
        physical = reference_parameters(F_bias = 0.0u"V/m")
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
        )
        hamiltonian = project_hamiltonians(problem, zeros(problem.numerical.N_z))
        shape = (
            problem.numerical.N_E,
            problem.numerical.N_k,
            problem.numerical.N_b,
            problem.numerical.N_b,
        )
        gamma = 0.02
        sigma_retarded = zeros(ComplexF64, shape)
        for e in axes(sigma_retarded, 1),
            m in axes(sigma_retarded, 2),
            a in axes(sigma_retarded, 3)

            sigma_retarded[e, m, a, a] = -0.5im * gamma
        end
        green_retarded, condition_number, dyson_scale = retarded_green(
            problem.grids.ε,
            hamiltonian,
            sigma_retarded;
            return_scale = true,
        )
        spectral = spectral_function(green_retarded)
        kBT = Float64(ustrip(u"eV", CODATA.kᴮₑᵥ * physical.Tᴸ)) / problem.scales.E₀_eV
        occupation = @. inv(exp(problem.grids.ε / kBT) + 1)
        sigma_lesser = zeros(ComplexF64, shape)
        for e in axes(sigma_lesser, 1),
            m in axes(sigma_lesser, 2),
            a in axes(sigma_lesser, 3)

            sigma_lesser[e, m, a, a] = im * occupation[e] * gamma
        end
        green_lesser = keldysh_green(green_retarded, sigma_lesser)
        green_greater = greater_green(green_retarded, green_lesser)
        sigma_greater = zeros(ComplexF64, shape)
        for e in axes(sigma_greater, 1),
            m in axes(sigma_greater, 2),
            a in axes(sigma_greater, 3)

            sigma_greater[e, m, a, a] = -im * (1 - occupation[e]) * gamma
        end
        green = GreenState(
            green_retarded,
            green_lesser,
            green_greater,
            spectral,
            condition_number,
            dyson_scale,
        )
        reservoir = SelfEnergyFamily(sigma_retarded, sigma_lesser, sigma_greater)

        dimensionless_flux = QCLNEGF._boundary_flux_complex(problem, reservoir, green)
        @test abs(dimensionless_flux) < 2e-13
        resolved = energy_resolved_current(problem, reservoir, green)
        @test maximum(abs, ustrip.(u"A/m^2/eV", resolved)) < 1e-3
        public_flux = boundary_flux(problem, reservoir, green)
        @test isfinite(ustrip(u"m^-2*s^-1", public_flux))
    end
end

end # independent suite
