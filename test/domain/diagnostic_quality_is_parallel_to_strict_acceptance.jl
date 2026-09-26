module Suite_T058
include("../support/common.jl")
include("../support/numerical_convergence.jl")

@testset "Diagnostic quality is parallel to strict acceptance" begin
    tolerances = SolverTolerances()
    disabled_policy =
        ConvergencePolicy(stagnation_window = 0, stagnation_relative_improvement = 0.0)
    disabled = scba_diagnostic_quality_assessment(
        _scba_iteration(10; r_K = 9e-4, r_Σ = 9e-4, r_λ = 9e-5),
        tolerances,
        disabled_policy,
    )
    @test !disabled.eligible
    @test !disabled.passed
    @test disabled.limiting_metric == :diagnostic_quality_disabled

    enabled_policy = ConvergencePolicy(
        stagnation_window = 0,
        stagnation_relative_improvement = 0.0,
        diagnostic_quality = DiagnosticQualityPolicy(enabled = true),
    )
    approximate_iteration = _scba_iteration(10; r_K = 9e-4, r_Σ = 9e-4, r_λ = 9e-5)
    @test !scba_convergence_assessment(approximate_iteration, tolerances, enabled_policy).passed
    approximate = scba_diagnostic_quality_assessment(
        approximate_iteration,
        tolerances,
        enabled_policy,
    )
    @test approximate.passed
    @test isempty(approximate.failed_metrics)
    state = QCLNEGF._scba_diagnostic_quality_state(2, approximate, enabled_policy)
    @test state.streak == 3
    @test state.quality == :approximate_fixed_point

    # Algebraic and causality gates stay strict; observables have an explicit approximate threshold.
    unsafe = scba_diagnostic_quality_assessment(
        _scba_iteration(
            10;
            r_K = 9e-4,
            r_Σ = 9e-4,
            r_λ = 9e-5,
            r_caus = 2e-10,
            r_Jchange = 2e-3,
        ),
        tolerances,
        enabled_policy,
    )
    @test !unsafe.passed
    @test :r_caus in unsafe.failed_metrics
    @test :r_Jchange in unsafe.failed_metrics
end

end # independent suite
