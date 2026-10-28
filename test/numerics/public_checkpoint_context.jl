module PublicCheckpointContextTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")
const N=QCLNEGF.QCLNumerics
struct StopContext <: Exception end

function caught(call)
    try
        call()
    catch error
        return error
    end
    return nothing
end
error_text(error)=error isa Exception ? sprint(showerror,error) : "no error"
numerical_payload(result)=(result.Uᴴ,result.n,result.scba,result.outer_history,result.report,
    result.status,result.converged,result.observables[:warnings],result.observables[:inner_working_policy])
function without_physical_markers(result)
    # Compare every numerical/decision field. Sampling of markers is checked separately below.
    rows=[Tuple(name===:physical_markers ? nothing : getfield(row,name)
        for name in fieldnames(SCBAIteration)) for row in result.scba.history]
    scba=Tuple(name===:history ? rows : getfield(result.scba,name) for name in fieldnames(SCBAResult))
    return Base.setindex(numerical_payload(result),scba,3)
end

@testset "Explicit sink admission rejects ambiguity before restore or compute" begin
    f=public_state_fixture()
    phases=Ref(0)
    production=QCLNEGF.with_production_options(f.production;phase_request=request->begin
        phases[]+=1
        request.task_width
    end)
    old=state->nothing
    context=(problem,state)->nothing
    for entry in (solve_production,solve_adaptive_production)
        dual=caught(()->entry(f.problem;options=f.options,production_options=production,
            execution=f.execution,checkpoint_sink=old,checkpoint_context_sink=context))
        @test dual isa ArgumentError
        @test occursin("mutually exclusive",error_text(dual))
        # A bad handle must not reach restore/hash validation before the keyword gate.
        bad_restart=caught(()->entry(f.problem;restart=:unvalidated,
            checkpoint_sink=old,checkpoint_context_sink=context))
        @test bad_restart isa ArgumentError
        @test occursin("mutually exclusive",error_text(bad_restart))
        missing_identity=caught(()->entry(f.problem;options=f.options,
            production_options=production,checkpoint_context_sink=context))
        @test missing_identity isa ArgumentError
        @test occursin("execution identity",error_text(missing_identity))
    end
    missing_kernels=caught(()->solve_adaptive_production(f.problem;options=f.options,
        production_options=production,execution=f.execution,
        domain_adaptation=DomainAdaptationPolicy(mode=:expand_energy_window),
        checkpoint_context_sink=context))
    @test missing_kernels isa ArgumentError
    @test occursin("kernel_options",error_text(missing_kernels))
    @test phases[]==0
end

@testset "Fixed and domain-none contexts preserve trajectory and detached ownership" begin
    problem=public_state_problem()
    # Explicit synthetic branch controls exercise several mixing steps, not scientific accuracy.
    tolerances=SolverTolerances(; (name=>1e100 for name in fieldnames(SolverTolerances))...)
    options=SolverOptions(α_Σ=0.1,max_scba=3,max_poisson=2,tolerances=tolerances,
        convergence=ConvergencePolicy(minimum_scba_iterations=3,required_consecutive_scba_passes=1,
            required_consecutive_poisson_passes=1,stagnation_window=0,
            stagnation_relative_improvement=0.0))
    production=ProductionOptions(worker_count=1,checkpoint_every_scba=1,checkpoint_every_outer=1,
        progress_every_scba=0,progress_every_outer=0)
    execution=public_state_execution(production)
    bind=(;options,production_options=production,execution,
        domain_adaptation=DomainAdaptationPolicy(),kernel_options=nothing)
    disabled=solve_production(problem;options,production_options=production)
    old_states=ExactSolverState[]
    legacy=solve_production(problem;options,production_options=production,execution,
        checkpoint_sink=state->push!(old_states,state))
    @test !isempty(old_states)
    states=ExactSolverState[]
    changed_snapshot=Ref(false)
    final_problem_id=Ref{UInt}(0)
    result=solve_production(problem;options,production_options=production,execution,
        checkpoint_context_sink=(current,state)->begin
            @test current===problem
            @test validate_exact_state(current,state;bind...)===nothing
            state.boundary===:final && (final_problem_id[]=objectid(current))
            # Mutate one retained snapshot while the actual solve is still active.
            if !changed_snapshot[] && state.boundary===:scba_accepted
                state.scba.green.Gᴿ[1]=987.0+2im
                state.continuation===nothing || state.continuation.candidate===nothing ||
                    (state.continuation.candidate.plus.Σᴿ[1]=654.0+3im)
                changed_snapshot[]=true
            end
            push!(states,state)
        end)
    @test changed_snapshot[]
    @test final_problem_id[]==objectid(result.problem)
    @test N._exact_digest(numerical_payload(result))==N._exact_digest(numerical_payload(legacy))
    @test N._exact_digest(without_physical_markers(result))==
        N._exact_digest(without_physical_markers(disabled))
    # Existing 5b policy: checkpoint_due forces terminal marker sampling at intermediate ν=2.
    # A source-isolated old-baseline probe reproduced exactly this sole no-sink difference.
    @test disabled.scba.history[2].physical_markers===nothing
    @test result.scba.history[2].physical_markers isa SCBAPhysicalMarkers
    for index in (1,3)
        @test N._exact_digest(result.scba.history[index].physical_markers)==
            N._exact_digest(disabled.scba.history[index].physical_markers)
    end
    trace(rows)=[(s.boundary,s.progress,s.identity,s.status) for s in rows]
    @test trace(states)==trace(old_states)
    @test length(result.scba.history)>=3
    @test Set(s.boundary for s in states)==Set((:scba_accepted,:outer_accepted,:final))
    original=result.scba.green.Gᴿ[1]
    last(states).scba.green.Gᴿ[1]=777.0+4im
    @test result.scba.green.Gᴿ[1]==original

    handle=restore_exact_state(problem,first(old_states);bind...)
    for entry in (solve_production,solve_adaptive_production)
        resumed_states=ExactSolverState[]
        resumed=entry(problem;restart=handle,checkpoint_context_sink=(current,state)->begin
            @test current===problem
            @test validate_exact_state(current,state;bind...)===nothing
            push!(resumed_states,state)
        end)
        @test !isempty(resumed_states)
        @test N._exact_digest(numerical_payload(resumed))==N._exact_digest(numerical_payload(result))
        @test last(resumed_states).boundary===:final
        @test_throws ArgumentError entry(problem;restart=handle,initial_Uᴴ=zeros(problem.numerical.N_z),
            checkpoint_context_sink=(current,state)->nothing)
    end
    none_states=ExactSolverState[]
    none_result=solve_adaptive_production(problem;options,production_options=production,execution,
        domain_adaptation=DomainAdaptationPolicy(),checkpoint_context_sink=(current,state)->begin
            @test current===problem
            push!(none_states,state)
        end)
    @test trace(none_states)==trace(old_states)
    @test N._exact_digest(numerical_payload(none_result))==N._exact_digest(numerical_payload(result))

    seed_before=N._exact_digest((result.Uᴴ,result.scba))
    warm=ExactSolverState[]
    solve_production(problem;options,production_options=production,execution,
        initial_Uᴴ=result.Uᴴ,initial_scba=result.scba,
        checkpoint_context_sink=(current,state)->begin
            @test current===problem
            push!(warm,state)
        end)
    @test !isempty(warm) && first(warm).progress.accepted_scba==1
    @test N._exact_digest((result.Uᴴ,result.scba))==seed_before

    @test_throws StopContext solve_production(problem;options,production_options=production,execution,
        checkpoint_context_sink=(current,state)->throw(StopContext()))
    # No arity guessing/fallback to the state-only legacy protocol.
    @test_throws MethodError solve_production(problem;options,production_options=production,execution,
        checkpoint_context_sink=state->throw(StopContext()))
end

@testset "Adaptive accepted states carry the actual rebuilt problem through pending decisions" begin
    f=public_state_fixture()
    numerical=NumericalParameters(; (name=>(name===:k_max ? 0.1u"nm^-1" :
        name===:E_min ? -1.0u"eV" : name===:E_max ? 1.0u"eV" :
        getfield(f.problem.numerical,name)) for name in fieldnames(NumericalParameters))...)
    problem=build_problem(;physical=reference_parameters(f_ion=1e-8),numerical,
        scales=f.problem.scales,scattering=f.problem.scattering,physical_models=f.problem.models)
    options=SolverOptions(max_scba=1,max_poisson=1)
    production=f.production
    execution=f.execution
    policy=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=1,
        maximum_energy_nodes=13,tail_threshold=1e-100)
    kernels=ProductionKernelOptions()
    bind=(;options,production_options=production,execution,domain_adaptation=policy,kernel_options=kernels)
    seen=Tuple{Int,Symbol,Symbol}[]
    cold=Ref{Union{Nothing,ExactSolverState}}(nothing)
    pending=Ref{Union{Nothing,ExactSolverState}}(nothing)
    final=Ref{Union{Nothing,ExactSolverState}}(nothing)
    final_problem_id=Ref{UInt}(0)
    result=solve_adaptive_production(problem;options,production_options=production,execution,
        domain_adaptation=policy,kernel_options=kernels,checkpoint_context_sink=(current,state)->begin
            @test validate_exact_state(current,state;bind...)===nothing
            @test current.numerical.N_E==size(state.scba.green.Gᴿ,1)
            push!(seen,(current.numerical.N_E,state.boundary,state.adaptation.next_action))
            if current.numerical.N_E!=problem.numerical.N_E
                @test current!==problem
                @test_throws ArgumentError validate_exact_state(problem,state;bind...)
                state.boundary===:domain_final && state.adaptation.next_action===:assess_measured_domain &&
                    (pending[]=state)
            end
            state.boundary===:domain_final && state.adaptation.next_action===:cold_rebuild && (cold[]=state)
            if state.boundary===:final
                final[]=state;final_problem_id[]=objectid(current)
            end
        end)
    @test Set(first(row) for row in seen)==Set((9,13))
    @test cold[]!==nothing && pending[]!==nothing && final[]!==nothing
    @test result.problem.numerical.N_E==13
    @test final_problem_id[]==objectid(result.problem)
    @test last(seen)==(13,:final,:finished)
    if cold[]!==nothing && pending[]!==nothing
        cold_handle=restore_exact_state(problem,cold[];bind...)
        resumed_domains=Int[]
        cold_result=solve_adaptive_production(problem;restart=cold_handle,
            checkpoint_context_sink=(current,state)->begin
                push!(resumed_domains,current.numerical.N_E)
                @test current.numerical.N_E==13
                @test validate_exact_state(current,state;bind...)===nothing
                @test_throws ArgumentError validate_exact_state(problem,state;bind...)
            end)
        @test !isempty(resumed_domains) && all(==(13),resumed_domains)
        @test N._exact_digest(numerical_payload(cold_result))==N._exact_digest(numerical_payload(result))
        pending_handle=restore_exact_state(result.problem,pending[];bind...)
        history=ScientificHistorySnapshot[]
        pending_states=ExactSolverState[]
        pending_result=solve_adaptive_production(result.problem;restart=pending_handle,
            history_observer=row->push!(history,row),checkpoint_context_sink=(current,state)->begin
                @test current===result.problem
                @test validate_exact_state(current,state;bind...)===nothing
                push!(pending_states,state)
            end)
        @test isempty(history) # Pending domain assessment does not run SCBA/Poisson again.
        @test length(pending_states)==1 && only(pending_states).boundary===:final
        @test pending_result.problem===result.problem
        # Legacy and context ports have identical complete payloads on the same pending branch.
        legacy_pending_states=ExactSolverState[]
        legacy_pending=solve_adaptive_production(result.problem;restart=pending_handle,
            checkpoint_sink=state->push!(legacy_pending_states,state))
        @test N._exact_digest(numerical_payload(pending_result))==
            N._exact_digest(numerical_payload(legacy_pending))
        @test [(s.boundary,s.progress,s.identity,s.status) for s in pending_states]==
            [(s.boundary,s.progress,s.identity,s.status) for s in legacy_pending_states]
        # Existing 5b terminal reconstruction retains policy evaluations but does not repeat
        # fresh-result descriptors. Check the precise difference, not a broad payload exclusion.
        @test N._exact_digest(Base.front(numerical_payload(pending_result)))==
            N._exact_digest(Base.front(numerical_payload(result)))
        original_policy=result.observables[:inner_working_policy]
        pending_policy=pending_result.observables[:inner_working_policy]
        @test Set(keys(original_policy))==Set(("mode","evaluations_since_resume",
            "strict_final_required","forcing_rule"))
        @test Set(keys(pending_policy))==Set(("evaluations_since_resume",))
        @test N._exact_digest(pending_policy["evaluations_since_resume"])==
            N._exact_digest(original_policy["evaluations_since_resume"])
    end
end
end
