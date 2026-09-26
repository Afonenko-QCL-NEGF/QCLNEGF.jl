module InexactInnerRequiresStrictFinal
include("../support/common.jl")
row(i, u, n) = OuterIteration(i, 0.0, u, n, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, [1.0], [1.0])

@testset "Raw outer residual controls inexact inner tolerance, never final gates" begin
    BN = QCLNEGF
    options = SolverOptions(
        convergence = ConvergencePolicy(
            mode = :adaptive_working,
            diagnostic_quality = DiagnosticQualityPolicy(enabled = true),
        ),
    )
    early = BN.inner_working_policy(options, OuterIteration[])
    @test early.stage === :inexact_inner
    @test early.options.tolerances === options.tolerances
    history = [row(1, 0.1, Inf), row(2, 0.01, 0.01), row(3, 1e-4, 2e-4)]
    later = BN.inner_working_policy(options, history)
    @test later.options.convergence.diagnostic_quality.keldysh_threshold ≈ 2e-5
    @test later.options.convergence.diagnostic_quality.self_energy_threshold ≈ 2e-5
    @test later.options.tolerances === options.tolerances
    regressed = BN.inner_working_policy(options, [history; row(4, 0.1, 0.1)])
    @test regressed.forcing_threshold == later.forcing_threshold
    strict = BN.inner_working_policy(
        options,
        [
            history;
            row(4, options.tolerances.r_U, options.tolerances.r_n)
        ],
    )
    @test strict.stage === :strict_final
    @test !strict.options.convergence.diagnostic_quality.enabled
    @test strict.options.tolerances === options.tolerances
    @test strict.options.max_scba == options.max_scba
    @test strict.options.max_poisson == options.max_poisson
    research = SolverOptions(
        convergence = ConvergencePolicy(
            mode = :research_continue,
            diagnostic_quality = DiagnosticQualityPolicy(enabled = true),
        ),
    )
    @test BN.inner_working_policy(research, history).options === research
end
end
