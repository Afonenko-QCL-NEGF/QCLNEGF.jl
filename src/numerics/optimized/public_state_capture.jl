mutable struct _ExactCaptureContext
    identity::ExactRestartIdentity
    options::SolverOptions
    production::ProductionOptions
    execution::ExactExecutionIdentity
    adaptation::Union{Nothing,DomainAdaptationState}
    sink::Function
    last_progress::Union{Nothing,SolverProgress}
    last_mixer_position::Union{Nothing,MixerPosition}
    domain_started_ns::Union{Nothing,UInt64}
    domain_elapsed_prior::Float64
end

"""Fixed-Hartree production solve. Exact continuation control is available only through a validated handle."""
function solve_scba_production(problem::NEGFProblem,Uᴴ::AbstractVector{<:Real};
    options::SolverOptions=baseline_options(),production_options::ProductionOptions=ProductionOptions(),
    cache::Union{Nothing,ProductionCache}=nothing,initial::Union{Nothing,SCBAResult}=nothing,
    consume_initial::Bool=false,iteration_callback::Union{Nothing,Function}=nothing,
    history_observer::Union{Nothing,Function}=nothing,resume_iterations::Bool=false,
    restart_options::SolverOptions=options)
    return _solve_scba_production(problem,Uᴴ;options,production_options,cache,initial,
        consume_initial,iteration_callback,history_observer,resume_iterations,restart_options)
end
function _exact_context(problem,options,production,execution,policy,kernel_options,sink;
    adaptation=nothing)
    identity=exact_restart_identity(problem,options,production,execution;
        domain_adaptation=policy,kernel_options)
    return _ExactCaptureContext(identity,options,production,execution,adaptation,sink,nothing,nothing,nothing,0.0)
end
function _adaptation_at_checkpoint(ctx)
    a=ctx.adaptation
    a===nothing && return nothing
    seconds=ctx.domain_started_ns===nothing ? a.elapsed_domain_seconds :
        ctx.domain_elapsed_prior+(time_ns()-ctx.domain_started_ns)/1e9
    values=Tuple(name===:elapsed_domain_seconds ? seconds : getfield(a,name)
        for name in fieldnames(DomainAdaptationState))
    # Borrow the immutable template until ExactSolverState performs its one owning copy.
    return DomainAdaptationState(values...)
end
function _acceptance_context(production;observations=nothing)
    raw=observations===nothing ? Dict{String,Any}() : get(observations,:acceptance_measurement_context,Dict())
    return AcceptanceContext(backend=Symbol(get(raw,"backend","production")),
        algorithms=production.algorithms,
        hilbert_roundoff_measured=get(raw,"hilbert_roundoff_measured",nothing),
        map_kernels_shared_with_solver=get(raw,"map_kernels_shared_with_solver",true),
        final_state_alignment=Symbol(get(raw,"final_state_alignment","producer_aligned")))
end
_records(rows)=Tuple(CoreRecord(row) for row in rows)
function _capture_state(ctx,problem,boundary,U,density,scba,outer,warnings,evaluations,progress,
    continuation,status;final_evidence=nothing,context=_acceptance_context(ctx.production))
    _validate_execution_runtime(ctx.execution)
    adaptation=_adaptation_at_checkpoint(ctx)
    # Conservative copy + canonical byte-buffer headroom; no host probe or per-event traversal.
    payload=(U,density,scba,outer,warnings,evaluations,continuation,final_evidence,adaptation)
    burst=Base.Checked.checked_add(Base.Checked.checked_mul(Base.summarysize(payload),4),32768)
    _production_phase_scope(
        token->begin
            state=ExactSolverState(;identity=ctx.identity,boundary,published_U=U,published_n=density,
                scba,mixer_position=ctx.last_mixer_position,
                family_order=Tuple((problem.kernels.enabled...,:embedding_plus,:embedding_minus)),
                outer_history=outer,progress,continuation,status,warnings=_records(warnings),
                policy_evaluations=_records(evaluations),adaptation,
                measurement_context=context,final_evidence)
            ctx.last_progress=progress
            ctx.sink(state)
            nothing
        end,
        ()->_production_phase_begin(ctx.production,:exact_state_snapshot,progress.accepted_scba;
            stage=:checkpoint,task_width=1,elastic=false,memory_burst_bytes=burst),
        (token,status)->_production_phase_end(ctx.production,token;status),
    )
    return nothing
end
function _capture_inner(ctx,problem,U,scba,outer,warnings,evaluations,μ,outer_streak,
    continuation,strict,diagnostic,position)
    ctx.last_mixer_position=position
    progress=SolverProgress(length(outer),μ,length(scba.history),outer_streak,strict,diagnostic)
    return _capture_state(ctx,problem,:scba_accepted,U,
        _electron_density_bar(problem,scba.green.Gˡ),scba,outer,warnings,evaluations,progress,
        continuation,:running_scba)
end
function _history_streaks(scba,options)
    strict=_trailing_pass_count(scba.history,row->scba_convergence_assessment(
        row,options.tolerances,options.convergence).passed)
    diagnostic=options.convergence.mode===:adaptive_working &&
        options.convergence.diagnostic_quality.enabled ?
        _trailing_pass_count(scba.history,row->scba_diagnostic_quality_assessment(
            row,options.tolerances,options.convergence).passed) : 0
    return strict,diagnostic
end
function _capture_outer(ctx,problem,U,next_U,density,scba,outer,warnings,evaluations,μ,streak,effective)
    strict,diagnostic=_history_streaks(scba,effective)
    done=streak>=ctx.options.convergence.required_consecutive_poisson_passes || μ>=ctx.options.max_poisson
    continuation=OuterContinuation(next_U,done ? :finalize : :next_outer)
    return _capture_state(ctx,problem,:outer_accepted,U,density,scba,outer,warnings,evaluations,
        SolverProgress(length(outer),μ,length(scba.history),streak,strict,diagnostic),continuation,:running)
end
const _EXACT_FINAL_OBSERVATIONS = (:acceptance_measurement_context,:final_fixed_point_audit,
    :charge_constraint,:restart_contract,:inner_working_policy,:domain_adaptation,
    :adaptation_checkpoint,:quality,:termination_reason,:physical_models,:effective_seed,
    :representation_coverage,:electron_electron,:lo_population)
function _final_evidence(result)
    observations=Dict{String,Any}(String(key)=>result.observables[key] for key in
        _EXACT_FINAL_OBSERVATIONS if haskey(result.observables,key))
    return FinalEvidenceSnapshot(report=result.report,original_converged=result.converged,
        original_status=result.status,observations=CoreRecord(observations))
end
function _capture_final(ctx,result,streak,evaluations;progress_override=nothing)
    outer=result.outer_history
    μ=isempty(evaluations) ? length(outer) : Int(get(last(evaluations),"outer_iteration",length(outer)))
    completed_before=μ<=length(outer) ? outer[1:max(0,μ-1)] : outer
    effective=inner_working_policy(ctx.options,completed_before).options
    strict,diagnostic=_history_streaks(result.scba,effective)
    boundary=ctx.adaptation===nothing || ctx.adaptation.next_action===:finished ? :final : :domain_final
    progress=progress_override===nothing ?
        SolverProgress(length(outer),μ,length(result.scba.history),streak,strict,diagnostic) : progress_override
    return _capture_state(ctx,result.problem,boundary,result.Uᴴ,result.n,result.scba,outer,
        get(result.observables,:warnings,Dict{String,Any}[]),evaluations,
        progress,
        nothing,result.status;final_evidence=_final_evidence(result),
        context=_acceptance_context(ctx.production;observations=result.observables))
end

function _public_observers(production,observation_sink,history_observer;exact=false)
    sequence=Ref(0); dropped=Ref(0)
    next_sequence()=(sequence[]+=1)
    old_sink=production.event_sink
    event_sink=observation_sink===nothing ? old_sink : function(event)
        owned=old_sink===nothing ? true : old_sink(event)
        outer=event.outer_iteration
        outer===nothing && production.outer_iteration>0 && (outer=production.outer_iteration)
        for metric in event.metrics
            metric.name===:outer_iteration && metric.value isa Integer &&
                (outer=metric.value>0 ? Int(metric.value) : nothing)
        end
        record=SolverEvent(event.action,event.stage,event.label,event.status,event.iteration,
            event.total,(event.metrics...,SolverMetric(:dropped_observations,dropped[],"1")),
            event.message,next_sequence(),outer,event.inner_iteration)
        try
            observation_sink(record)
        catch error
            error isa InterruptException && rethrow()
            dropped[]+=1
        end
        owned
    end
    old_phase=production.phase_request
    phase=exact && old_phase!==nothing ? function(request)
        width=old_phase(request)
        width==request.task_width || throw(ArgumentError(
            "execution.phase_request changed numerical width under exact identity"))
        width
    end : old_phase
    production=with_production_options(production;event_sink,phase_request=phase)
    scientific=history_observer===nothing ? nothing :
        function(problem,domain_revision)
            # Bind the represented domain once, including history-only invocations.
            g,b=problem.grids,problem.basis
            discretization_id=_representation_digest(g.x,g.wˣ,g.ε,g.wᴱ,g.κ,g.wᵏ,g.qᶻ,g.φ)
            basis_id=_representation_digest(b.Φ,b.H₀,b.M⁻¹)
            return function(stage,μ,row,decision,strict,diagnostic,action)
                snapshot=ScientificHistorySnapshot(;action,stage,domain_revision,discretization_id,basis_id,
                    outer_iteration=μ,
                    inner_iteration=stage===:scba ? row.ν : nothing,sequence=next_sequence(),
                    converged=decision.converged,status=decision.status,quality=decision.quality,
                    strict_streak=strict,diagnostic_streak=diagnostic,row)
                history_observer(snapshot)
                nothing
            end
        end
    return production,scientific,dropped
end

function _restore_binding(problem,restart)
    _state_require(restart isa RestoredSolverState,"restart","requires a validated Core handle")
    validate_exact_state(problem,restart.state;options=restart.options,
        production_options=restart.production_options,execution=restart.execution,
        domain_adaptation=restart.domain_adaptation,kernel_options=restart.kernel_options)
    return deepcopy(restart.state)
end
function _restart_numerical_keywords(state)
    c=state.continuation
    U=c isa OuterContinuation ? c.next_U : state.published_U
    return (;initial_Uᴴ=U,initial_published_U=state.published_U,initial_scba=state.scba,
        initial_outer_history=state.outer_history,initial_warnings=Dict{String,Any}[_core_plain(w) for w in state.warnings],
        initial_policy_evaluations=Dict{String,Any}[_core_plain(row) for row in state.policy_evaluations],
        initial_outer_streak=state.progress.outer_strict_streak,
        saved_continuation=c isa SCBAContinuation ? c : nothing,
        saved_mixer_position=state.mixer_position,
        resume_scba=c isa SCBAContinuation,consume_initial=true)
end
function _public_solve_production(problem::NEGFProblem;restart=nothing,options=nothing,
    production_options=nothing,execution=nothing,checkpoint_sink=nothing,
    observation_sink=nothing,history_observer=nothing,kwargs...)
    state=nothing
    if restart!==nothing
        isempty(kwargs) && options===nothing && production_options===nothing && execution===nothing ||
            throw(ArgumentError("restart does not accept overriding inputs or legacy initial_*"))
        state=_restore_binding(problem,restart)
        _state_require(state.boundary in (:scba_accepted,:outer_accepted) && state.adaptation===nothing,
            "restart.boundary","requires fixed-domain continuation")
        options=restart.options; production_options=restart.production_options; execution=restart.execution
    else
        options===nothing && (options=baseline_options())
        production_options===nothing && (production_options=ProductionOptions())
    end
    exact=state!==nothing || checkpoint_sink!==nothing
    exact && state===nothing && get(kwargs,:resume_scba,false) &&
        throw(ArgumentError("exact continuation requires a validated handle with mixer position"))
    exact && execution===nothing && throw(ArgumentError("execution identity is required for exact snapshots"))
    production,scientific,dropped=_public_observers(production_options,observation_sink,history_observer;exact)
    ctx=checkpoint_sink===nothing ? nothing : _exact_context(problem,options,production,
        execution,DomainAdaptationPolicy(),nothing,checkpoint_sink)
    ctx===nothing || state===nothing || (ctx.last_mixer_position=state.mixer_position)
    initial=state===nothing ? (;kwargs...) : _restart_numerical_keywords(state)
    result=_solve_production(problem;options,production_options=production,exact_capture=ctx,
        scientific_sink=scientific===nothing ? nothing : scientific(problem,0),initial...)
    result.observables[:observation_dropped]=dropped[]
    return result
end
