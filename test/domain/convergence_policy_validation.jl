module Suite_T059
include("../support/common.jl")
include("../support/numerical_convergence.jl")

@testset "Convergence policy validation" begin
    valid = SolverOptions()
    @test QCLNEGF._check_options(valid) === valid
    @test_throws ArgumentError QCLNEGF._check_options(
        SolverOptions(
            convergence = ConvergencePolicy(required_consecutive_scba_passes = 0),
        ),
    )
    @test_throws ArgumentError QCLNEGF._check_options(
        SolverOptions(convergence = ConvergencePolicy(stagnation_window = 3)),
    )
    @test_throws ArgumentError QCLNEGF._check_options(
        SolverOptions(
            convergence = ConvergencePolicy(
                stagnation_window = 0,
                stagnation_relative_improvement = 0.1,
            ),
        ),
    )
    @test_throws ArgumentError QCLNEGF._check_options(
        SolverOptions(
            convergence = ConvergencePolicy(
                diagnostic_quality = DiagnosticQualityPolicy(keldysh_threshold = 0.0),
            ),
        ),
    )
    @test_throws ArgumentError QCLNEGF._check_options(
        SolverOptions(
            convergence = ConvergencePolicy(
                diagnostic_quality = DiagnosticQualityPolicy(
                    required_consecutive_passes = 2,
                ),
            ),
        ),
    )
end

end # independent suite
