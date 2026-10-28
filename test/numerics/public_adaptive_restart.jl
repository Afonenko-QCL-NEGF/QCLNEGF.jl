module PublicAdaptiveRestartTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")
struct StopPendingDomain <: Exception end

@testset "Pending measured-domain resume retains cost recorded before its checkpoint" begin
    problem=public_state_problem()
    policy=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=0,maximum_energy_nodes=31)
    options=SolverOptions(α_Σ=0.1,max_scba=1,max_poisson=1)
    initializations=Ref(0)
    production=ProductionOptions(worker_count=1,checkpoint_every_scba=0,checkpoint_every_outer=0,
        progress_every_scba=0,progress_every_outer=0,phase_request=request->begin
            if request.name===:initialization
                initializations[]+=1
                sleep(2.0) # Explicit bounded engineering delay inside the measured original domain.
            end
            request.task_width
        end)
    execution=public_state_execution(production)
    kernels=ProductionKernelOptions()
    saved=Ref{Union{Nothing,ExactSolverState}}(nothing)
    @test_throws StopPendingDomain solve_adaptive_production(problem;options,
        production_options=production,execution,domain_adaptation=policy,kernel_options=kernels,
        checkpoint_sink=state->begin
            if state.boundary===:domain_final
                saved[]=state
                throw(StopPendingDomain())
            end
        end)
    @test saved[]!==nothing
    if saved[]!==nothing
        handle=restore_exact_state(problem,saved[];options,production_options=production,
            execution,domain_adaptation=policy,kernel_options=kernels)
        resume_started=time_ns()
        result=solve_adaptive_production(problem;restart=handle)
        resume_wall=(time_ns()-resume_started)*1e-9
        @test initializations[]==1
        seconds=only(result.observables[:domain_adaptation]["revisions"])["seconds"]
        println("Recorded domain cost: original initialization delay=2s; saved=",
            saved[].adaptation.elapsed_domain_seconds,"; resume_wall=",resume_wall,"; seconds=",seconds)
        @test seconds>=2.0
        @test saved[].adaptation.elapsed_domain_seconds>=2.0
        @test seconds>=saved[].adaptation.elapsed_domain_seconds
        # The new domain interval is contained in this external invocation interval;
        # a second charge of the saved prior must not be hidden in the revision cost.
        @test seconds<=saved[].adaptation.elapsed_domain_seconds+resume_wall+1e-6
    end
end

@testset "Domain final resumes its pending adaptive decision rather than its solve" begin
    f=public_state_fixture()
    policy=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=0,
        maximum_energy_nodes=31)
    kernels=ProductionKernelOptions()
    identity=exact_restart_identity(f.problem,f.options,f.production,f.execution;
        domain_adaptation=policy,kernel_options=kernels)
    adaptation=DomainAdaptationState(policy=policy,kernel_options=kernels,expansions=0,elapsed_domain_seconds=0.0,
        current_problem=identity.problem,current_problem_recipe=ProblemRebuildRecipe(
            physical=f.problem.physical,numerical=f.problem.numerical,scales=f.problem.scales,
            scattering=f.problem.scattering,models=f.problem.models,
            algorithms=f.production.algorithms,kernel_options=kernels),
        revisions=(),next_action=:assess_measured_domain,
        decision_id=0,decision_committed=false,next_problem_recipe=nothing)
    evidence=FinalEvidenceSnapshot(report=ConvergenceReport(false,Dict(:r_K=>NaN),["postmortem domain"]),
        original_converged=false,original_status=:scba_invalid_candidate,
        observations=CoreRecord(Dict("quality"=>"unconverged")))
    state=ExactSolverState(identity=identity,boundary=:domain_final,
        published_U=f.state.published_U,published_n=f.state.published_n,scba=public_state_failed_scba(f.state.scba),
        mixer_position=f.state.mixer_position,family_order=f.state.family_order,outer_history=f.state.outer_history,
        progress=f.state.progress,continuation=nothing,status=:scba_invalid_candidate,
        warnings=f.state.warnings,policy_evaluations=(),adaptation=adaptation,
        measurement_context=f.state.measurement_context,final_evidence=evidence)
    handle=restore_exact_state(f.problem,state;options=f.options,production_options=f.production,
        execution=f.execution,domain_adaptation=policy,kernel_options=kernels)
    history=ScientificHistorySnapshot[]
    states=ExactSolverState[]
    result=try
        solve_adaptive_production(f.problem;restart=handle,
            checkpoint_sink=row->push!(states,row),history_observer=row->push!(history,row))
    catch error
        error
    end
    @test typeof(result)===NEGFSolution
    if result isa NEGFSolution
        @test isempty(history) # Completed domain did not emit another SCBA/Poisson row.
        @test result.Uᴴ == state.published_U
        @test result.scba.green.Gᴿ == state.scba.green.Gᴿ
        @test length(result.observables[:domain_adaptation]["revisions"])==1
        @test result.observables[:domain_adaptation]["expansions"]==0
        @test last(states).boundary===:final
        @test last(states).adaptation.next_action===:finished
        @test last(states).adaptation.decision_committed
        @test last(states).adaptation.kernel_options==kernels
        @test hasproperty(last(states).adaptation,:current_problem_recipe)
        if hasproperty(last(states).adaptation,:current_problem_recipe)
            @test last(states).adaptation.current_problem_recipe.numerical===f.problem.numerical
        end
    end
    different=ProductionKernelOptions(relative_tolerance=nextfloat(kernels.relative_tolerance))
    @test_throws ArgumentError restore_exact_state(f.problem,state;options=f.options,
        production_options=f.production,execution=f.execution,domain_adaptation=policy,
        kernel_options=different)
end

@testset "Cold rebuild consumes the bound recipe without replacing kernel policy" begin
    f=public_state_fixture()
    kernels=ProductionKernelOptions(relative_tolerance=1e-6,initial_lookup_nodes=17,
        maximum_lookup_nodes=33,impurity_angular_tolerance=2e-7)
    recipe=ProblemRebuildRecipe(physical=f.problem.physical,numerical=f.problem.numerical,
        scales=f.problem.scales,scattering=f.problem.scattering,models=f.problem.models,
        algorithms=f.production.algorithms,kernel_options=kernels)
    seen=Ref{Any}(nothing)
    builder=(;kwargs...)->begin
        seen[]=(;kwargs...)
        (;problem=f.problem)
    end
    result=QCLNEGF.QCLNumerics._rebuild_from_recipe(recipe;builder)
    @test result===f.problem
    @test seen[].kernel_options===kernels
    @test seen[].algorithms===f.production.algorithms
    @test seen[].physical===recipe.physical
    @test seen[].numerical===recipe.numerical
    @test_throws MethodError solve_scba_production(f.problem,zeros(f.problem.numerical.N_z);
        saved_continuation=:unvalidated)
end
end
