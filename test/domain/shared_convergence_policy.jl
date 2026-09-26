module Suite_T056
include("../support/common.jl")
include("../support/numerical_convergence.jl")

@testset "Shared convergence policy" begin
    tolerances = SolverTolerances()
    policy = ConvergencePolicy(
        minimum_scba_iterations = 2,
        stagnation_window = 0,
        stagnation_relative_improvement = 0.0,
    )

    too_early = scba_convergence_assessment(_scba_iteration(1), tolerances, policy)
    @test !too_early.eligible
    @test !too_early.passed
    @test isempty(too_early.failed_metrics)

    accepted = scba_convergence_assessment(_scba_iteration(2), tolerances, policy)
    @test accepted.eligible
    @test accepted.passed
    # r_D, r_Jchange, and r_population all have ratio 0.1 here.  The shared
    # policy resolves an exact tie by retaining the first declared gate.
    @test accepted.limiting_metric == :r_D
    @test accepted.limiting_ratio ≈ 0.1

    rejected = scba_convergence_assessment(
        _scba_iteration(3; r_K = 2e-8, r_Σ = 3e-8),
        tolerances,
        policy,
    )
    @test !rejected.passed
    @test rejected.limiting_metric == :r_Σ
    @test rejected.limiting_ratio ≈ 3.0
    @test rejected.failed_metrics == [:r_K, :r_Σ]

    nonfinite =
        scba_convergence_assessment(_scba_iteration(4; r_K = NaN), tolerances, policy)
    @test !nonfinite.passed
    @test nonfinite.limiting_metric == :r_K
    @test isinf(nonfinite.limiting_ratio)

    outer = OuterIteration(
        2,
        1e-12,
        1e-8,
        1e-8,
        1e-11,
        1e-7,
        1e-5,
        1e-5,
        1e-12,
        1.0,
        [1.0],
        [1.0],
    )
    @test poisson_convergence_assessment(outer, tolerances, policy).passed
end

end # independent suite
