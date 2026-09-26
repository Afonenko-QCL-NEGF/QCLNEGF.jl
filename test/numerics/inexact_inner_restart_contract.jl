module InexactInnerRestartContract
include("../support/common.jl")

@testset "Adaptive inner restart retains the immutable enclosing policy" begin
    BN=QCLNEGF
    problem=build_problem(
        numerical = tutorial_numerics(),
        physical = reference_parameters(Tᴸ = 70u"K"),
        scattering = ScatteringOptions(
            LO = false,
            acoustic = true,
            impurity = false,
            IFR = false,
            alloy = false,
        ),
    )
    base=SolverOptions(
        α_Σ = 0.2,
        max_scba = 3,
        convergence = ConvergencePolicy(
            mode = :adaptive_working,
            stagnation_window = 0,
            stagnation_relative_improvement = 0.0,
            diagnostic_quality = DiagnosticQualityPolicy(enabled = true),
        ),
    )
    outer=[OuterIteration(2, 0.0, 1e-4, 2e-4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, [1.0], [1.0])]
    effective=BN.inner_working_policy(base, outer).options
    production=ProductionOptions(worker_count = 1, checkpoint_every_scba = 1)
    saved=Ref{Union{Nothing,SCBAResult}}(nothing)
    uninterrupted=solve_scba_production(
        problem,
        zeros(problem.numerical.N_z);
        options = effective,
        restart_options = base,
        production_options = production,
        iteration_callback = state ->
            isempty(state.history) || last(state.history).ν!=1 ? nothing :
            (saved[]=deepcopy(state)),
    )
    @test saved[] !== nothing
    @test saved[].restart_contract==BN._solver_restart_contract(base, production.algorithms)
    @test_throws ArgumentError solve_scba_production(
        problem,
        zeros(problem.numerical.N_z);
        options = effective,
        production_options = production,
        initial = saved[],
        resume_iterations = true,
    )
    resumed=solve_scba_production(
        problem,
        zeros(problem.numerical.N_z);
        options = effective,
        restart_options = base,
        production_options = production,
        initial = saved[],
        resume_iterations = true,
    )
    @test resumed.status==uninterrupted.status
    @test resumed.green.Gᴿ ≈ uninterrupted.green.Gᴿ rtol=2e-13
    @test resumed.green.Gˡ ≈ uninterrupted.green.Gˡ rtol=2e-13
    @test [row.r_Σ for row in resumed.history] ≈ [row.r_Σ for row in uninterrupted.history] rtol=2e-13
    partial=NEGFSolution(
        problem,
        base,
        zeros(problem.numerical.N_z),
        zeros(problem.numerical.N_z),
        resumed,
        OuterIteration[],
        Dict{Symbol,Any}(),
        ConvergenceReport(false, Dict(:r_sum=>0.2, :r_U=>Inf), String[]),
        false,
        :max_scba,
    )
    assessment=solution_scientific_assessment(partial)
    @test !assessment["iterative_converged"]
    @test !assessment["scientific_accepted"]
    @test assessment["discretization_verified"]=="not_measured"
    @test assessment["metrics"]["r_sum"]["status"]=="fail"
    @test assessment["metrics"]["r_U"]["status"]=="error"
    @test assessment["metrics"]["r_K"]["status"]=="not_measured"
    @test assessment["metrics"]["r_K"]["value"]===nothing
end
end
