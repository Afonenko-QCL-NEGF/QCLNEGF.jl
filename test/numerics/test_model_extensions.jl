module ModelExtensionsTests
using Test
using Unitful
using LinearAlgebra
using QCLNEGF

let M = QCLNEGF.PhysicalModelExtensions
    @testset "Screening dimension and thermodynamic limits" begin
        mass = 0.067 * CODATA.m₀
        ε = 12.9 * CODATA.ε₀
        debye = M.DebyeScreening3D(1e20u"m^-3", 70u"K", 12.9)
        κ = M.screening_wavenumber(debye)
        # Debye length follows independently by linearizing Poisson-Boltzmann.
        λ = sqrt(ε * CODATA.kᴮ * 70u"K" / (CODATA.e^2 * 1e20u"m^-3"))
        @test ustrip(Unitful.NoUnits, κ * λ) ≈ 1
        @test M.screening_wavenumber(M.DebyeScreening3D(4e20u"m^-3", 70u"K", 12.9)) ≈ 2κ
        @test M.screening_wavenumber(M.DebyeScreening3D(1e20u"m^-3", 280u"K", 12.9)) ≈ κ / 2
        classical = M.ThomasFermiScreening2D(1e8u"m^-2", 70u"K", mass, 2, 12.9)
        classical_limit = CODATA.e^2 * 1e8u"m^-2" / (2ε * CODATA.kᴮ * 70u"K")
        @test M.screening_wavenumber(classical) ≈ classical_limit rtol = 1e-5
        degenerate = M.ThomasFermiScreening2D(1e20u"m^-2", 70u"K", mass, 2, 12.9)
        @test M.screening_wavenumber(degenerate) ≈ CODATA.e^2 * mass / (2π * ε * CODATA.ħ^2)
        @test iszero(
            M.screening_wavenumber(
                M.ThomasFermiScreening2D(0u"m^-2", 70u"K", mass, 2, 12.9),
            ),
        )
        fixed2 = M.FixedScreening(1e8u"m^-1", 12.9, 2)
        fixed3 = M.FixedScreening(1e8u"m^-1", 12.9, 3)
        @test M.screened_coulomb_potential(fixed2, 0u"m^-1") ≈
              CODATA.e^2 / (2ε * 1e8u"m^-1")
        @test M.screened_coulomb_potential(fixed3, 0u"m^-1") ≈
              CODATA.e^2 / (ε * 1e16u"m^-2")
        @test M.model_identity(fixed2) != M.model_identity(fixed3)
        @test_throws ArgumentError M.FixedScreening(1u"m^-1", 12.9, 1)
        @test_throws ArgumentError M.DebyeScreening3D(1u"m^-3", 0u"K", 12.9)
        @test_throws Unitful.DimensionError M.DebyeScreening3D(1u"m^-2", 70u"K", 12.9)
        @test_throws DomainError M.screened_coulomb_potential(
            M.FixedScreening(0u"m^-1", 12.9, 2),
            0u"m^-1",
        )
    end

    @testset "Kane dispersion, state counting and consistent velocity" begin
        mass = 0.067 * CODATA.m₀
        parabolic = M.KaneDispersion(mass, 0u"eV^-1", 2)
        kane = M.KaneDispersion(mass, 10u"eV^-1", 2) # manufactured, not a GaAs fit
        k = 1e8u"m^-1"
        e0 = uconvert(u"eV", CODATA.ħ^2 * k^2 / (2mass))
        e = M.kane_energy(kane, k)
        @test M.kane_energy(parabolic, k) ≈ e0
        @test e * (1 + kane.nonparabolicity * e) ≈ e0
        @test e < e0
        # Plane-wave counting is independent of the energy parameterization.
        @test M.kane_state_count_2d(kane, e) ≈ 2k^2 / (4π)
        h = 1e-4 * k
        numerical_velocity =
            (M.kane_energy(kane, k + h) - M.kane_energy(kane, k - h)) / (2h * CODATA.ħ)
        @test M.kane_velocity(kane, k) ≈ numerical_velocity rtol = 1e-7
        edges = [0, 0.01, 0.05, 0.1] .* u"eV"
        weights = M.kane_energy_bin_weights_2d(kane, edges)
        @test sum(weights) ≈ M.kane_state_count_2d(kane, edges[end])
        for j in eachindex(weights)
            @test weights[j] ≈
                  M.kane_dos_2d(kane, (edges[j] + edges[j+1]) / 2) * (edges[j+1] - edges[j])
        end
        rotation = [cos(0.4) -sin(0.4); sin(0.4) cos(0.4)]
        kinetic = rotation * Diagonal([0.01, 0.07]) * rotation' .* u"eV"
        kinetic = (kinetic + kinetic') / 2
        transformed = ustrip.(u"eV", M.kane_kinetic_operator(kane, kinetic))
        @test transformed + 10transformed^2 ≈ ustrip.(u"eV", kinetic)
        @test M.kane_kinetic_operator(parabolic, kinetic) ≈ kinetic
        @test_throws ArgumentError M.KaneDispersion(mass, -1u"eV^-1", 2)
        @test_throws DomainError M.kane_kinetic_operator(kane, [-1.0 0; 0 1] .* u"eV")
        @test_throws ArgumentError M.kane_energy_bin_weights_2d(kane, [0, 0] .* u"eV")
    end

    @testset "LO kinetics: detailed balance and energy accounting" begin
        bath = M.HotLOPhononBath(36.7u"meV", 70u"K", 5u"ps")
        n0 = M.equilibrium_lo_population(bath)
        exponent = ustrip(Unitful.NoUnits, bath.phonon_energy / (CODATA.kᴮₑᵥ * 70u"K"))
        factors = M.lo_scattering_factors(n0)
        @test factors.absorption / factors.emission ≈ exp(-exponent)
        @test M.lo_population_temperature(bath, n0) ≈ 70u"K"
        b = 1e12u"s^-1"
        a = b * exp(-exponent)
        @test M.steady_lo_population(bath, a, b).occupation ≈ n0.occupation
        hot = M.steady_lo_population(bath, 1e11u"s^-1", 2e11u"s^-1")
        @test hot.occupation > n0.occupation
        balance = M.lo_energy_balance(bath, hot, 1e11u"s^-1", 2e11u"s^-1")
        @test balance.electron_to_lo ≈ balance.lo_to_bath
        @test abs(ustrip(u"W", balance.stored_energy_rate)) < 1e-24
        initial = M.LOPopulation(2)
        cooled = M.advance_lo_population(bath, initial, 0u"s^-1", 0u"s^-1", 10u"ps")
        @test cooled.occupation ≈ n0.occupation + (2 - n0.occupation) * exp(-2)
        first = M.advance_lo_population(bath, initial, a, b, 1u"ps")
        second = M.advance_lo_population(bath, first, a, b, 2u"ps")
        combined = M.advance_lo_population(bath, initial, a, b, 3u"ps")
        @test second.occupation ≈ combined.occupation
        @test_throws DomainError M.steady_lo_population(bath, 1e12u"s^-1", 0u"s^-1")
        @test_throws ArgumentError M.LOPopulation(-0.01)
        @test_throws ArgumentError M.HotLOPhononBath(36.7u"meV", 70u"K", 0u"ps")
    end

    @testset "Pauli binary collisions conserve physical moments" begin
        energies = [0.01, 0.07, 0.03, 0.05] .* u"eV"
        wavevectors = [(1e8, 0), (-1e8, 0), (0, 1e8), (0, -1e8)]
        wavevectors = [v .* u"m^-1" for v in wavevectors]
        channel = M.BinaryCollisionChannel((1, 2, 3, 4), 1e12u"s^-1")
        model = M.ElectronCollisionModel(
            energies,
            wavevectors,
            [channel];
            rate_provenance = "manufactured reversible channel",
        )
        f = [0.8, 0.7, 0.2, 0.1]
        df = M.electron_collision_rhs(model, f)
        @test df[1] < 0u"s^-1" && df[3] > 0u"s^-1"
        # Project moments directly, without using the tested helper as oracle.
        @test sum(df) == 0u"s^-1"
        @test abs(ustrip(u"eV/s", sum(energies .* df))) < 1e-4
        @test sum(v[1] * d for (v, d) in zip(wavevectors, df)) == 0u"m^-1/s"
        @test sum(v[2] * d for (v, d) in zip(wavevectors, df)) == 0u"m^-1/s"
        # Fermi detailed balance follows from conserved E1+E2=E3+E4.
        fermi = [
            1 / (
                1 + exp(ustrip(Unitful.NoUnits, (e - 0.04u"eV") / (CODATA.kᴮₑᵥ * 70u"K")))
            ) for e in energies
        ]
        @test maximum(abs, ustrip.(u"s^-1", M.electron_collision_rhs(model, fermi))) < 1e-3
        @test all(iszero, M.electron_collision_rhs(model, ones(4)))
        @test all(iszero, M.electron_collision_rhs(model, zeros(4)))
        boundary = M.electron_collision_rhs(model, [1, 1, 0, 0])
        @test all(ustrip.(u"s^-1", boundary[1:2]) .<= 0)
        @test all(ustrip.(u"s^-1", boundary[3:4]) .>= 0)
        @test M.electron_collision_moments(model, df).number_rate == 0u"s^-1"
        bad_energies = copy(energies)
        bad_energies[1] += 1u"meV"
        @test_throws ArgumentError M.ElectronCollisionModel(
            bad_energies,
            wavevectors,
            [channel];
            rate_provenance = "fixture",
        )
        reverse = M.BinaryCollisionChannel((3, 4, 1, 2), 1e12u"s^-1")
        @test_throws ArgumentError M.ElectronCollisionModel(
            energies,
            wavevectors,
            [channel, reverse];
            rate_provenance = "fixture",
        )
        @test_throws ArgumentError M.electron_collision_rhs(model, [1.01, 0, 0, 0])
        @test_throws ArgumentError M.BinaryCollisionChannel((1, 1, 2, 3), 1u"s^-1")
    end

    @testset "Optical power sign and peak-field normalization" begin
        ω = 2π * 3e12u"s^-1"
        field = 2u"V/m"
        χ = 0.1 + 0.02im
        power = M.optical_power_density(χ, ω, field)
        conductivity = CODATA.ε₀ * ω * imag(χ)
        # Joule heating with RMS field = peak/sqrt(2).
        @test power ≈ conductivity * (field / sqrt(2))^2
        @test power > 0u"W/m^3"
        @test M.optical_power_density(conj(χ), ω, field) == -power
        @test M.optical_power_density(χ, ω, 2field) == 4power
        @test iszero(M.optical_power_density(0.1, ω, field))
        @test_throws ArgumentError M.optical_power_density(χ, -ω, field)
    end
end

end # independent suite
