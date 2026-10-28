function _rebuild_from_recipe(recipe::ProblemRebuildRecipe;builder=build_configured_scattering_problem)
    built=builder(physical=recipe.physical,numerical=recipe.numerical,scales=recipe.scales,
        scattering=recipe.scattering,physical_models=recipe.models,algorithms=recipe.algorithms,
        kernel_options=recipe.kernel_options)
    return built.problem
end
function _domain_revisions(rows)
    return Tuple(map(rows) do row
        source=Symbol(row["source"])
        T=source===:preflight_estimate ? DomainPreflightEvidence :
          source===:measured_state ? DomainMeasuredEvidence :
          throw(ArgumentError("adaptation.revision.source is unknown"))
        expected=Set(("revision","source",String.(fieldnames(T))...))
        Set(keys(row))==expected || throw(ArgumentError("adaptation.revision: missing or extra evidence fields"))
        evidence=T(; (name=>row[String(name)] for name in fieldnames(T))...)
        DomainRevision(;revision=Int(row["revision"]),source,evidence)
    end)
end
function _domain_revision_plain(row::DomainRevision)
    evidence=row.evidence
    record=Dict{String,Any}(String(name)=>deepcopy(getfield(evidence,name))
        for name in fieldnames(typeof(evidence)))
    record["revision"]=row.revision; record["source"]=String(row.source)
    return record
end
function _problem_recipe(problem,algorithms,kernels)
    return ProblemRebuildRecipe(physical=problem.physical,numerical=problem.numerical,
        scales=problem.scales,scattering=problem.scattering,models=problem.models,
        algorithms=algorithms,kernel_options=kernels)
end
function _adaptation_snapshot(policy,kernels,identity,current,algorithms,revisions,expansions,action,recipe;
    elapsed_domain_seconds)
    return DomainAdaptationState(;policy=policy,kernel_options=kernels,expansions,elapsed_domain_seconds,
        current_problem=identity.problem,current_problem_recipe=_problem_recipe(current,algorithms,kernels),
        revisions=_domain_revisions(revisions),
        next_action=action,decision_id=length(revisions),
        decision_committed=action in (:cold_rebuild,:finished),next_problem_recipe=recipe)
end

function _terminal_state_solution(problem,state,options,production)
    observations=_collect_observables(problem,state.scba,state.published_U)
    for pair in state.final_evidence.observations.entries
        observations[Symbol(first(pair))]=_core_plain(last(pair))
    end
    observations[:warnings]=Dict{String,Any}[_core_plain(row) for row in state.warnings]
    observations[:restart_contract]=_solver_restart_contract(options,production.algorithms)
    c=state.measurement_context
    observations[:acceptance_measurement_context]=Dict{String,Any}(
        "backend"=>String(c.backend),"hilbert_roundoff_measured"=>c.hilbert_roundoff_measured,
        "map_kernels_shared_with_solver"=>c.map_kernels_shared_with_solver,
        "final_state_alignment"=>String(c.final_state_alignment),"fresh_candidate_required"=>true)
    observations[:inner_working_policy]=Dict{String,Any}("evaluations_since_resume"=>
        Dict{String,Any}[_core_plain(row) for row in state.policy_evaluations])
    return NEGFSolution(problem,options,state.published_U,state.published_n,state.scba,
        state.outer_history,observations,deepcopy(state.final_evidence.report),
        state.final_evidence.original_converged,state.status)
end

"""Continue a declared adaptive execution, including a completed domain's pending decision."""
function solve_adaptive_production(problem::NEGFProblem;restart=nothing,options=nothing,
    production_options=nothing,execution=nothing,checkpoint_sink=nothing,
    checkpoint_context_sink::Union{Nothing,Function}=nothing,
    observation_sink=nothing,history_observer=nothing,domain_adaptation=nothing,
    kernel_options=nothing,kwargs...)
    sink=_checkpoint_callback(checkpoint_sink,checkpoint_context_sink)
    state=nothing
    if restart!==nothing
        isempty(kwargs) && options===nothing && production_options===nothing && execution===nothing &&
            domain_adaptation===nothing && kernel_options===nothing ||
            throw(ArgumentError("adaptive restart does not accept overriding inputs or legacy initial_*"))
        state=_restore_binding(problem,restart)
        _state_require(state.boundary in (:scba_accepted,:outer_accepted,:domain_final),
            "restart.boundary","adaptive execution is already global final")
        options=restart.options; production_options=restart.production_options; execution=restart.execution
        domain_adaptation=restart.domain_adaptation; kernel_options=restart.kernel_options
    else
        options===nothing && (options=reference_production_solver_options())
        production_options===nothing && (production_options=ProductionOptions())
        domain_adaptation===nothing && (domain_adaptation=DomainAdaptationPolicy())
    end
    if domain_adaptation.mode===:none
        state===nothing || return _public_solve_production(problem;restart,checkpoint_context_sink=sink,
            observation_sink,history_observer)
        return _public_solve_production(problem;options,production_options,execution,
            checkpoint_context_sink=sink,observation_sink,history_observer,kwargs...)
    end
    exact=state!==nothing || sink!==nothing
    exact && state===nothing && get(kwargs,:resume_scba,false) &&
        throw(ArgumentError("exact continuation requires a validated handle with mixer position"))
    exact && execution===nothing && throw(ArgumentError("execution identity is required for adaptive snapshots"))
    exact && kernel_options===nothing && throw(ArgumentError("kernel_options is required for exact rebuild"))
    kernel_options===nothing && (kernel_options=ProductionKernelOptions())
    production,scientific,dropped=_public_observers(production_options,observation_sink,history_observer;exact)
    contexts=Ref{Union{Nothing,_ExactCaptureContext}}(nothing)
    factory=sink===nothing ? nothing : function(current,revisions,expansions,started,elapsed_prior)
        ctx=_exact_context(current,options,production,execution,domain_adaptation,kernel_options,sink)
        ctx.domain_started_ns=started
        ctx.domain_elapsed_prior=elapsed_prior
        ctx.adaptation=_adaptation_snapshot(domain_adaptation,kernel_options,ctx.identity,current,production.algorithms,
            revisions,expansions,:assess_measured_domain,nothing;elapsed_domain_seconds=elapsed_prior)
        if state!==nothing && state.identity.problem==ctx.identity.problem
            ctx.last_progress=state.progress
            ctx.last_mixer_position=state.mixer_position
        end
        contexts[]=ctx
        return ctx
    end
    transition=sink===nothing ? nothing : function(latest,current,revisions,expansions,action,recipe)
        ctx=contexts[]
        ctx===nothing && throw(ArgumentError("adaptive transition has no bound capture context"))
        ctx.adaptation=_adaptation_snapshot(domain_adaptation,kernel_options,ctx.identity,current,production.algorithms,
            revisions,expansions,action,recipe;elapsed_domain_seconds=ctx.domain_elapsed_prior)
        evaluations=get(get(latest.observables,:inner_working_policy,Dict()),
            "evaluations_since_resume",Dict{String,Any}[])
        progress=ctx.last_progress
        progress===nothing && throw(ArgumentError("adaptive transition has no accepted progress"))
        _capture_final(ctx,latest,progress.outer_strict_streak,evaluations;progress_override=progress)
        nothing
    end
    initial=(;kwargs...)
    current=problem
    if state!==nothing
        adaptation=state.adaptation
        _state_require(adaptation!==nothing,"state.adaptation","missing")
        revisions=Dict{String,Any}[_domain_revision_plain(row) for row in adaptation.revisions]
        expansions=adaptation.expansions
        if state.boundary===:domain_final && adaptation.next_action===:cold_rebuild
            current=_rebuild_from_recipe(adaptation.next_problem_recipe)
            expansions+=1
            initial=(;initial_adaptation=Dict("revisions"=>revisions,"expansions"=>expansions))
        elseif state.boundary===:domain_final
            _state_require(adaptation.next_action===:assess_measured_domain,
                "state.adaptation.next_action","unsupported domain-final continuation")
            latest=_terminal_state_solution(problem,state,options,production)
            initial=(;initial_domain_solution=latest,
                initial_elapsed_domain_seconds=adaptation.elapsed_domain_seconds,
                initial_adaptation=Dict("revisions"=>revisions,"expansions"=>expansions))
        else
            initial=(;_restart_numerical_keywords(state)...,
                initial_elapsed_domain_seconds=adaptation.elapsed_domain_seconds,
                initial_adaptation=Dict("revisions"=>revisions,"expansions"=>expansions))
        end
    end
    result=_solve_adaptive_production(current;options,production_options=production,domain_adaptation,
        kernel_options,scientific_sink_factory=scientific,exact_capture_factory=factory,
        adaptive_transition_sink=transition,initial...)
    result.observables[:observation_dropped]=dropped[]
    return result
end
