module PairedOccupationNormalization
include("../support/common.jl")
const BN = QCLNEGF
const RO = BN.QCLReferenceOperators

@testset "Fixed-number map preserves the matrix Pauli cone without clipping" begin
    # A nearly full small block reproduces unilateral-scaling failure. The
    # second block supplies empty capacity. No eigenvalues are truncated.
    Gn = reshape(ComplexF64[0.999999, 0.1], 2, 1, 1, 1)
    Gp = 1 .- Gn
    target = 1.11
    legacy = Gn .* (target/sum(real, Gn))
    @test minimum(real.(1 .- legacy)) < 0
    lesser, greater = im .* Gn, -im .* Gp
    lambda =
        RO._normalize_keldysh_pair!(lesser, greater, sum(real, Gn), sum(real, Gp), target)
    @test sum(real, -im .* lesser) ≈ target rtol=4eps(Float64)
    @test minimum(real.(-im .* lesser)) >= 0
    @test minimum(real.(im .* greater)) >= 0
    @test im .* (greater-lesser) ≈ Gn+Gp atol=4eps(Float64)
    @test lambda < 1
    @test abs(lambda-1) > 1e-8 # enforcing charge does not certify a raw fixed point

    # Excess number transfers occupied weight to empty weight; fixed N is identity.
    for requested in (0.7, 1.099999)
        lesser, greater = im .* Gn, -im .* Gp
        RO._normalize_keldysh_pair!(
            lesser,
            greater,
            sum(real, Gn),
            sum(real, Gp),
            requested,
        )
        @test sum(real, -im .* lesser) ≈ requested rtol=4eps(Float64)
        @test minimum(real.(im .* greater)) >= 0
        @test im .* (greater-lesser) ≈ Gn+Gp atol=4eps(Float64)
    end
    @test_throws DomainError RO._paired_number_coefficients(1.0, 0.1, 1.2)
    @test_throws DomainError RO._paired_number_coefficients(1.0, -0.1, 0.8)
end

@testset "Actual production Dyson / Keldysh update keeps direct empty-state weight" begin
    scattering = ScatteringOptions(
        LO = false,
        acoustic = false,
        impurity = false,
        IFR = false,
        alloy = false,
    )
    problem = build_problem(numerical = tutorial_numerics(), scattering = scattering)
    h = project_hamiltonians(problem, zeros(problem.numerical.N_z))
    n = problem.numerical
    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    sigmaR, sigmaL, sigmaG = (zeros(ComplexF64, shape) for _ = 1:3)
    gamma = 0.02
    for e = 1:n.N_E, m = 1:n.N_k, a = 1:n.N_b
        sigmaR[e, m, a, a] = -0.5im*gamma
    end
    GR, _, _ = BN._retarded_green_production(problem.grids.ε, h, sigmaR)
    A = spectral_function(GR)
    capacity = number_functional(im .* A, problem.grids, problem.physical.g_s)
    target = BN._scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    raw_number = 0.99target
    critical = zeros(ComplexF64, shape)
    critical[1, 1, :, :] .= im .* A[1, 1, :, :]
    critical_capacity = number_functional(critical, problem.grids, problem.physical.g_s)
    f_full = 0.999999
    f_other = (raw_number-f_full*critical_capacity)/(capacity-critical_capacity)
    @test 0 < f_other < 1
    for e = 1:n.N_E, m = 1:n.N_k, a = 1:n.N_b
        f = e == 1 && m == 1 ? f_full : f_other
        sigmaL[e, m, a, a] = im*f*gamma
        sigmaG[e, m, a, a] = -im*(1-f)*gamma
    end
    embedding = SelfEnergyFamily(sigmaR, sigmaL, sigmaG)
    families = Dict{Symbol,SelfEnergyFamily}()
    corrected, lambda, total =
        BN._current_green_production(problem, h, families, embedding, ProductionOptions())
    legacy, legacy_lambda, _ = BN._current_green_production(
        problem,
        h,
        families,
        embedding,
        ProductionOptions(
            algorithms = AlgorithmOptions(occupation_normalization = :scalar_lesser),
        ),
    )
    @test lambda ≈ 0.99 rtol=2e-13
    @test legacy_lambda == lambda
    @test number_functional(corrected.Gˡ, problem.grids, problem.physical.g_s) ≈ target rtol=2e-13
    @test eigmin(Hermitian(im .* legacy.Gᵍ[1, 1, :, :])) < -1e-6
    @test eigmin(Hermitian(im .* corrected.Gᵍ[1, 1, :, :])) >= 0
    @test corrected.A ≈ im .* (corrected.Gᵍ-corrected.Gˡ) rtol=1e-13
    @test BN._green_positivity(corrected, total) < 1e-12
    @test BN._keldysh_residual(corrected, total) > 1e-6

    serial =
        BN.cavity_embedding_self_energy(problem, total, h; periods = 2, worker_count = 1)
    parallel =
        BN.cavity_embedding_self_energy(problem, total, h; periods = 2, worker_count = 2)
    for index in eachindex(serial), component in (:Σᴿ, :Σˡ, :Σᵍ)
        @test getfield(serial[index], component) == getfield(parallel[index], component)
    end

    # A raw negative out-scattering component remains visible after the map.
    invalid = BN._copy_family(total)
    invalid.Σᵍ[1, 1, 1, 1] = 0.01im
    metrics = production_residual_suite(problem.grids.ε, h, corrected, invalid, invalid)
    @test metrics.r_PSD > 0.1
end

@testset "Anderson guards the entire self-energy Keldysh triple" begin
    family = SelfEnergyFamily(
        fill(-0.5im, 1, 1, 1, 1),
        fill(0.4im, 1, 1, 1, 1),
        fill(-0.6im, 1, 1, 1, 1),
    )
    @test BN._anderson_causal([family])
    # Retarded remains causal and identity remains exact; incoming rate is negative.
    bad = SelfEnergyFamily(
        fill(-0.5im, 1, 1, 1, 1),
        fill(-0.1im, 1, 1, 1, 1),
        fill(-1.1im, 1, 1, 1, 1),
    )
    @test !BN._anderson_causal([bad])
    @test_throws DomainError BN._anderson_mix!(
        BN._AndersonWorkspace(),
        [family],
        [bad],
        AlgorithmOptions(mixing = :anderson),
        1.0,
        ProductionOptions(),
    )
end

@testset "Prescribed SPPA bath does not require zero electronic power" begin
    t = SolverTolerances()
    metrics =
        Dict(:r_energy_electron_electron=>0.5, :ee_energy_conservation_applicable=>0.0)
    limits = BN._solution_validation_limits(t, metrics)
    @test !haskey(limits, :r_energy_electron_electron)
    @test haskey(limits, :r_power)
    metrics[:ee_energy_conservation_applicable] = 1.0
    @test BN._solution_validation_limits(t, metrics)[:r_energy_electron_electron] ==
          t.r_power
end

@testset "SCBA retains first/worst/last matrix payloads with all scalar witnesses" begin
    history = SCBAIteration[]
    retention = BN.QCLNumerics._SCBAWitnessRetention(history)
    ratios = [1.0, 0.5, 0.2, 2.0, 0.1, 0.05]
    for (index, ratio) in enumerate(ratios)
        witness = BN.SCBAPhysicsWitness(
            :raw_unoccupied,
            1,
            1,
            -ratio,
            1.0,
            0.0,
            ratio,
            ratio,
            ratio,
            0.0,
            Dict("Gp_raw_direct_dimensionless"=>fill(ComplexF64(-ratio), 1, 1)),
            ComplexF64[1],
        )
        base = SCBAIteration(
            index,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            ratio,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            NaN,
            NaN,
            NaN,
            NaN,
            nothing,
            nothing,
        )
        push!(
            history,
            SCBAIteration(
                (
                    name === :witness ? witness : getfield(base, name) for
                    name in fieldnames(SCBAIteration)
                )...,
            ),
        )
        BN.QCLNumerics._retain_scba_witness_payloads!(history, retention)
    end
    @test [row.witness.ratio for row in history] == ratios
    @test findall(row -> !isempty(row.witness.matrices), history) == [1, 4, 6]
    resumed = BN.QCLNumerics._SCBAWitnessRetention(history)
    @test (resumed.first, resumed.worst, resumed.last) == (1, 4, 6)
end
end
