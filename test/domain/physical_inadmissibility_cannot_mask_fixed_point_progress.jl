module Suite_T064
include("../support/common.jl")
include("../support/numerical_convergence.jl")

@testset "Physical inadmissibility cannot mask fixed-point progress" begin
    tolerances = SolverTolerances()
    policy = ConvergencePolicy(
        mode = :adaptive_working,
        minimum_scba_iterations = 1,
        stagnation_window = 4,
        stagnation_relative_improvement = 0.05,
    )
    # Regression for the all_legacy diagnostic: r_PSD is almost exactly one,
    # while raw Keldysh/self-energy residuals keep decreasing.
    improving = [
        _scba_iteration(i; r_Σ = (5-i)*1e-3, r_K = (5-i)*2e-3, r_PSD = 0.999999965) for
        i = 1:4
    ]
    @test all(
        row ->
            scba_convergence_assessment(row, tolerances, policy).limiting_metric === :r_PSD,
        improving,
    )
    @test !scba_convergence_stagnated(improving, tolerances, policy)
    @test QCLNEGF._scba_stopping_reason(improving, tolerances, policy) === nothing
    # A true fixed-point plateau remains a bounded failure, never acceptance.
    plateau =
        [_scba_iteration(i; r_Σ = 0.038, r_K = 0.0704, r_PSD = 0.999999965) for i = 1:4]
    @test scba_convergence_stagnated(plateau, tolerances, policy)
    @test QCLNEGF._scba_stopping_reason(plateau, tolerances, policy) === :stagnated
    # Once the fixed point passes throughout the window, persistent physical
    # invalidity has an explicit diagnosis rather than fictitious convergence.
    blocked = [_scba_iteration(i; r_PSD = 0.999999965) for i = 1:4]
    @test QCLNEGF._scba_stopping_reason(blocked, tolerances, policy) === :quality_blocked
    @test !scba_convergence_stagnated(blocked, tolerances, policy)
    @test !scba_convergence_assessment(last(blocked), tolerances, policy).passed
    options = SolverOptions(
        convergence = ConvergencePolicy(
            mode = :adaptive_working,
            minimum_scba_iterations = 1,
            stagnation_window = 4,
            stagnation_relative_improvement = 0.05,
            diagnostic_quality = DiagnosticQualityPolicy(enabled = true),
        ),
    )
    @test !QCLNEGF.scba_approximate_acceptance(blocked, options)
    recovering = [_scba_iteration(i; r_PSD = (5-i)*1e-3) for i = 1:4]
    @test QCLNEGF._scba_stopping_reason(recovering, tolerances, policy) === nothing
    nonfinite = copy(plateau)
    nonfinite[2]=_scba_iteration(2; r_K = Inf, r_PSD = 1.0)
    @test QCLNEGF._scba_stopping_reason(nonfinite, tolerances, policy) === nothing
end

end # independent suite
