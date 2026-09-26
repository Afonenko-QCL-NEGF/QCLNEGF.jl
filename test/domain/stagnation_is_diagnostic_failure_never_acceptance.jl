module Suite_T057
include("../support/common.jl")
include("../support/numerical_convergence.jl")

@testset "Stagnation is diagnostic failure, never acceptance" begin
    tolerances = SolverTolerances()
    policy =
        ConvergencePolicy(stagnation_window = 4, stagnation_relative_improvement = 0.05)
    stalled = [_scba_iteration(i; r_Σ = (2.00 - 0.01i) * 1e-8) for i = 1:4]
    @test scba_convergence_stagnated(stalled, tolerances, policy)
    @test !scba_convergence_assessment(last(stalled), tolerances, policy).passed

    improving = [_scba_iteration(i; r_Σ = (3.0 - 0.45i) * 1e-8) for i = 1:4]
    @test !scba_convergence_stagnated(improving, tolerances, policy)
    disabled =
        ConvergencePolicy(stagnation_window = 0, stagnation_relative_improvement = 0.0)
    @test !scba_convergence_stagnated(stalled, tolerances, disabled)
end

end # independent suite
