module PublicObservationEquivalenceTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")

function numerical_payload(result)
    # No timing, dropped-event counters or caller-selected telemetry enter this comparison.
    return (result.Uᴴ,result.n,result.scba,result.outer_history,result.status,result.converged,
        result.report,result.observables[:warnings],result.observables[:inner_working_policy])
end

@testset "Technical observers and diagnostic mode preserve the numerical trajectory" begin
    problem=public_state_problem()
    # Explicit branch-control data only: these thresholds are not scientific acceptance evidence.
    tolerances=SolverTolerances(; (field=>1e100 for field in fieldnames(SolverTolerances))...)
    options=SolverOptions(α_Σ=0.1,max_scba=3,max_poisson=2,tolerances=tolerances,
        convergence=ConvergencePolicy(minimum_scba_iterations=2,
            required_consecutive_scba_passes=2,required_consecutive_poisson_passes=1,
            stagnation_window=0,stagnation_relative_improvement=0.0))
    production=ProductionOptions(worker_count=1,progress_every_scba=0,progress_every_outer=0)
    disabled=solve_production(problem;options,production_options=production)
    modes=(:normal,:diagnostic,:failing)
    for mode in modes
        events=SolverEvent[]
        rows=ScientificHistorySnapshot[]
        attempts=Ref(0)
        collector=function(event)
            attempts[]+=1
            push!(events,event)
            mode===:failing && error("synthetic technical collector failure")
            return nothing
        end
        observed=solve_production(problem;options,
            production_options=QCLNEGF.QCLNumerics.with_production_options(production;phase_timing=mode===:diagnostic),
            observation_sink=collector,history_observer=row->push!(rows,row))
        @test QCLNEGF.QCLNumerics._exact_digest(numerical_payload(observed))==
            QCLNEGF.QCLNumerics._exact_digest(numerical_payload(disabled))
        @test length(observed.scba.history)==3
        @test length(observed.outer_history)==2
        @test !isempty(events) && !isempty(rows)
        coverage=representation_coverage(problem)
        @test all(row->row.domain_revision==0 && row.discretization_id==coverage.discretization_id &&
            row.basis_id==coverage.basis_id,rows)
        @test all(row->row.discretization_id===first(rows).discretization_id &&
            row.basis_id===first(rows).basis_id,rows) # One cached binding, no per-row digest.
        @test issorted([event.sequence for event in events])
        @test allunique([event.sequence for event in events])
        @test all(event->event.metrics isa Tuple,events)
        headers=Dict{Tuple,Tuple}()
        for snapshot in rows
            snapshot.stage===:scba || continue
            key=(snapshot.domain_revision,snapshot.discretization_id,snapshot.basis_id,
                snapshot.stage,snapshot.outer_iteration,snapshot.inner_iteration)
            header=(snapshot.converged,snapshot.status,snapshot.quality,
                snapshot.strict_streak,snapshot.diagnostic_streak)
            if snapshot.action===:append
                headers[key]=header
            else
                @test haskey(headers,key)
                @test header==headers[key]
            end
        end
        @test any(row->row.action===:witness_update,rows)
        scba=filter(event->event.stage===:scba && event.action in (:begin,:end,:progress),events)
        @test all(event->event.outer_iteration in (1,2),scba)
        @test any(event->event.action===:end && event.stage===:scba,events)
        @test observed.observables[:observation_dropped]==(mode===:failing ? attempts[] : 0)
        @test !solution_scientific_assessment(observed)["scientific_accepted"]
        lastrow=last(filter(row->row.stage===:outer && row.action===:append,rows))
        original=copy(observed.outer_history[end].density)
        lastrow.row.density[1]=1234.0
        @test observed.outer_history[end].density==original
    end
end

@testset "Observation wrapper preserves existing stage ownership callbacks" begin
    legacy=SolverEvent[]
    production=ProductionOptions(event_sink=event->(push!(legacy,event);true))
    wrapped,_,_=QCLNEGF.QCLNumerics._public_observers(production,_ -> nothing,nothing)
    owned=QCLNEGF.QCLNumerics._begin_solver_stage(wrapped,:scba;iteration=1)
    @test owned
    QCLNEGF.QCLNumerics._end_solver_stage(wrapped,:scba,owned)
    @test [event.action for event in legacy]==[:begin,:end]
end
end
