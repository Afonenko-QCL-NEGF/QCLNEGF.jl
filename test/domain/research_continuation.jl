module ResearchContinuationContracts
include("../support/common.jl")

row(i, r; lambda = r, psd = 0.3) = SCBAIteration(
    i,
    1e-14,
    1e-12,
    r,
    r,
    lambda,
    1+lambda,
    psd,
    1e-12,
    1e-14,
    1e-5,
    1e-5,
    1.0,
    NaN,
    NaN,
    NaN,
    NaN,
    nothing,
    nothing,
)
function trajectory(values; lambda = nothing)
    [
        row(i, value; lambda = lambda === nothing ? value : lambda[i]) for
        (i, value) in enumerate(values)
    ]
end
@testset "Research trend contract version 2 (synthetic, no empirical evidence implied)" begin
    policy = ConvergencePolicy(
        mode = :research_continue,
        trend = ResearchTrendPolicy(window = 32),
        diagnostic_quality = DiagnosticQualityPolicy(enabled = true),
    )
    options = SolverOptions(max_scba = 2000, convergence = policy)
    smooth = exp.(range(log(1.0), log(1e-5); length = 2000))
    plateau = vcat(exp.(range(log(1.0), log(1e-5); length = 500)), fill(1e-5, 1500))
    for values in (smooth, plateau, plateau .* (1 .+ 0.2sin.(1:2000)))
        history = trajectory(values)
        decision = research_continuation_decision(history, options)
        @test decision.action === :continue_poisson
        @test decision.selected_state === :last
        @test decision.keldysh.reduction > 0.9
        @test !QCLNEGF.scba_approximate_acceptance(history[1:1999], options)
        @test QCLNEGF._scba_stopping_reason(
            history[1:1999],
            options.tolerances,
            policy,
        ) === nothing
        @test research_continuation_decision(history[1:1999], options).action === :iterate
    end
    @test research_continuation_decision(trajectory(fill(1.0, 2000)), options).action ===
          :stop
    @test research_continuation_decision(trajectory(reverse(smooth)), options).action ===
          :stop
    large_oscillation = repeat([1e-5, 1.0], 1000)
    @test research_continuation_decision(trajectory(large_oscillation), options).action ===
          :stop
    nonfinite = copy(plateau)
    nonfinite[end]=NaN
    @test research_continuation_decision(trajectory(nonfinite), options).reason ===
          :unusable_state
    @test research_continuation_decision(trajectory(plateau), options; usable = false).reason ===
          :unusable_state
    # r_lambda must have its own trend, even when both other raw residuals fall.
    @test research_continuation_decision(
        trajectory(plateau; lambda = fill(0.5, 2000)),
        options,
    ).action === :stop
    # A transient crossing of lambda=1 is not a stable normalization fixed
    # point. A later full window within the useful band must not be rejected
    # merely for being larger than that accidental historical minimum.
    crossing = vcat(
        exp.(range(log(0.8), log(1e-7); length = 500)),
        fill(1e-7, 32),
        fill(1.5e-5, 1468),
    )
    crossing_history = trajectory(plateau; lambda = crossing)
    decision = research_continuation_decision(crossing_history, options)
    @test decision.normalization.regression > policy.trend.max_regression
    @test decision.action === :continue_poisson
    @test !scba_convergence_assessment(last(crossing_history), options.tolerances, policy).passed
    @test !QCLNEGF.scba_approximate_acceptance(crossing_history, options)
    # Even starting exactly inside the useful band need not create an error
    # just to demonstrate its reduction. This remains a research decision.
    @test QCLNEGF.residual_trend(vcat(zeros(32), fill(1.5e-5, 32)), policy.trend).promising
    # Recovery above the band still needs a bounded regression; an isolated
    # spike above the band cannot be hidden by a small window average.
    backsliding = copy(crossing)
    backsliding[(end-31):end] .= 3e-4
    @test research_continuation_decision(
        trajectory(plateau; lambda = backsliding),
        options,
    ).action === :stop
    spike = copy(crossing)
    spike[end] = 3e-4
    @test research_continuation_decision(trajectory(plateau; lambda = spike), options).action ===
          :stop
    @test research_continuation_decision(crossing_history[1:1999], options).action ===
          :iterate
    @test_throws ArgumentError ResearchTrendPolicy(1, 0.5, 1e-4, 2.0, 10.0)
    @test_throws ArgumentError ResearchTrendPolicy(window = 32, promising_threshold = NaN)
end
end
