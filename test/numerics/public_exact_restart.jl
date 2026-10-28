module PublicExactRestartTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")

@testset "Exact state validates before compute and preserves diagnostic bits" begin
    @test isdefined(QCLNEGF,:validate_exact_state)
    if isdefined(QCLNEGF,:validate_exact_state)
        f=public_state_fixture()
        bind=(options=f.options,production_options=f.production,execution=f.execution,
            domain_adaptation=f.policy,kernel_options=nothing)
        @test validate_exact_state(f.problem,f.state;bind...) === nothing
        restored=restore_exact_state(f.problem,f.state;bind...)
        @test restored isa RestoredSolverState
        @test isnan(f.state.scba.history[1].lambda_change)
        @test isinf(f.state.scba.history[1].r_Jchange)
        @test isnan(f.state.warnings[1]["unmeasured"].value)
        changed=deepcopy(f.problem)
        changed.grids.ε[1]=nextfloat(changed.grids.ε[1])
        @test_throws ArgumentError validate_exact_state(changed,f.state;bind...)
        bad=ExactExecutionIdentity(;
            (name=>(name===:julia ? "1.13.1" : getfield(f.execution,name))
             for name in fieldnames(ExactExecutionIdentity))...)
        @test_throws ArgumentError exact_restart_identity(f.problem,f.options,f.production,bad;
            domain_adaptation=f.policy,kernel_options=nothing)
        f.state.scba.green.Gᴿ[1]=2.0+3im
        @test_throws ArgumentError validate_exact_state(f.problem,f.state;bind...)
    end
end

@testset "Accepted SCBA and outer snapshots resume without repeating an iteration" begin
    problem=public_state_problem()
    # Declared synthetic branch-control fixture, following the existing packaged restart test.
    # No threshold/default of a scientific calculation is changed by this fixture.
    tolerances=SolverTolerances(; (field=>1e100 for field in fieldnames(SolverTolerances))...)
    options=SolverOptions(α_Σ=0.1,max_scba=4,max_poisson=3,tolerances=tolerances,
        convergence=ConvergencePolicy(minimum_scba_iterations=4,required_consecutive_scba_passes=1,
            required_consecutive_poisson_passes=1,stagnation_window=0,
            stagnation_relative_improvement=0.0))
    production=ProductionOptions(worker_count=1,checkpoint_every_scba=1,
        checkpoint_every_outer=1,progress_every_scba=0,progress_every_outer=0)
    execution=public_state_execution(production)
    states=ExactSolverState[]
    observed=ScientificHistorySnapshot[]
    result=try
        solve_production(problem;options,production_options=production,execution,
            checkpoint_sink=state->push!(states,state),history_observer=row->push!(observed,row))
    catch error
        error
    end
    @test typeof(result) === NEGFSolution
    if result isa NEGFSolution
        @test Set(state.boundary for state in states) == Set((:scba_accepted,:outer_accepted,:final))
        @test !isempty(observed)
        @test all(row->row isa ScientificHistorySnapshot,observed)
        targets=Set((row.domain_revision,row.discretization_id,row.basis_id,row.stage,
            row.outer_iteration,row.inner_iteration) for row in observed)
        for state in states
            state.boundary===:final && continue
            handle=restore_exact_state(problem,state;options,production_options=production,
                execution,domain_adaptation=DomainAdaptationPolicy(),kernel_options=nothing)
            rows=ScientificHistorySnapshot[]
            resumed=solve_production(problem;restart=handle,history_observer=row->push!(rows,row))
            @test all(row->(row.domain_revision,row.discretization_id,row.basis_id,row.stage,
                row.outer_iteration,row.inner_iteration) in targets,rows)
            @test resumed.status == result.status
            @test resumed.Uᴴ == result.Uᴴ
            @test resumed.scba.green.Gᴿ == result.scba.green.Gᴿ
            @test resumed.scba.green.Gˡ == result.scba.green.Gˡ
            for field in fieldnames(SCBAIteration)[1:17]
                @test isequal([getfield(r,field) for r in resumed.scba.history],
                    [getfield(r,field) for r in result.scba.history])
            end
            @test [r.μ for r in resumed.outer_history] == [r.μ for r in result.outer_history]
            @test all(row->!(row.action===:append && row.stage===:scba && row.outer_iteration==state.progress.current_outer &&
                row.inner_iteration<=state.progress.accepted_scba),rows)
        end
        final=only(filter(state->state.boundary===:final,states))
        handle=restore_exact_state(problem,final;options,production_options=production,
            execution,domain_adaptation=DomainAdaptationPolicy(),kernel_options=nothing)
        readonly=restore_completed_solution(handle;revalidation=:selected_map,
            independent_oracle=:if_supported)
        @test readonly.Uᴴ==result.Uᴴ
        @test readonly.n==result.n
        for field in fieldnames(GreenState)
            @test isequal(getfield(readonly.scba.green,field),getfield(result.scba.green,field))
        end
        for field in fieldnames(SCBAIteration)[1:17]
            @test isequal([getfield(r,field) for r in readonly.scba.history],
                [getfield(r,field) for r in result.scba.history])
        end
        @test length(readonly.outer_history)==length(result.outer_history)
        @test readonly.observables[:restoration_evidence].revalidation.status===:measured
        @test readonly.observables[:restoration_evidence].original.original_status===result.status
        @test_throws ArgumentError solve_production(problem;restart=handle)
    end
end
end
