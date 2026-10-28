function _state_require(condition,path,reason)
    condition || throw(ArgumentError("$path: $reason"))
end
function _state_array(array,shape,path;finite=true)
    size(array)==shape || throw(DimensionMismatch("$path: expected $shape, got $(size(array))"))
    finite && _state_require(all(isfinite,array),path,"non-finite mandatory array")
end
function _state_family(family,shape,path;finite=true)
    for name in (:Σᴿ,:Σˡ,:Σᵍ)
        _state_array(getfield(family,name),shape,"$path.$name";finite)
    end
end
_state_payload(state::ExactSolverState)=
    Tuple(getfield(state,name) for name in fieldnames(ExactSolverState)[1:end-1])

function _validate_mixer_position(position,mixer,accepted,algorithms)
    _state_require(position.applied_mixes==max(accepted-1,0),
        "state.mixer_position.applied_mixes","differs from accepted pre-mix boundary")
    _state_require(position.last_coefficient_fallback in (:none,:coefficient_solve,:coefficient_bound),
        "state.mixer_position.last_coefficient_fallback","unknown fallback")
    _state_require(position.applied_mixes>0 || position.last_coefficient_fallback===:none,
        "state.mixer_position.last_coefficient_fallback","fallback before the first mix")
    if mixer.method===:linear
        _state_require(isempty(position.retained_mix_indices) && position.last_reset_mix==0 &&
            position.last_reset_reason===:none && position.last_coefficient_fallback===:none,
            "state.mixer_position","linear mixer has Anderson receipt")
    else
        applied=position.applied_mixes
        reason=position.last_reset_reason
        _state_require((reason===:none && position.last_reset_mix==0) ||
            (reason in (:residual_increase,:causality_guard) && 1<=position.last_reset_mix<=applied),
            "state.mixer_position.last_reset_mix","unknown or inconsistent reset")
        first_possible=reason===:none ? 1 :
            position.last_reset_mix+(reason===:causality_guard ? 1 : 0)
        first_retained=max(first_possible,applied-algorithms.anderson_history_depth+1)
        expected=Tuple(first_retained:applied)
        _state_require(position.retained_mix_indices==expected,
            "state.mixer_position.retained_mix_indices","missing or inconsistent mix coverage")
    end
    _state_require(length(mixer.states)==length(mixer.residuals)==length(position.retained_mix_indices),
        "state.mixer.history","missing payload for declared mix coverage")
    return nothing
end

function _recipe_energy_window(grid,problem,path)
    lower,upper,nodes=grid
    _state_require(isfinite(lower) && isfinite(upper) && isfinite(upper-lower) &&
        upper-lower>2*_electronvolts(problem.numerical.M_E) && nodes>=3,
        path,"invalid recorded eV window or energy-node count")
end
function _recipe_energy_nodes(ratio,path)
    _state_require(isfinite(ratio) && ratio>0,path,"energy-step ratio is not finite and positive")
    try
        return ceil(Int,ratio)+1
    catch error
        error isa InexactError || error isa OverflowError || rethrow()
        throw(ArgumentError("$path: energy-node count is not representable"))
    end
end
function _scaled_energy_preimage(value,scale,reference)
    # Enclose x for the actual preflight conversion fl(fl(x*scale)+reference).
    # Each inverse operation encloses its IEEE rounding cell; no rtol or node slack.
    product_lower=prevfloat(prevfloat(value)-reference)
    product_upper=nextfloat(nextfloat(value)-reference)
    return (prevfloat(prevfloat(product_lower)/scale),
        nextfloat(nextfloat(product_upper)/scale))
end
function _validate_recipe_energy_transition(old,new,problem,path)
    _recipe_energy_window(old,problem,path)
    _recipe_energy_window(new,problem,path)
    lower,upper,nodes=old
    next_lower,next_upper,next_nodes=new
    scale,reference=problem.scales.E₀_eV,_electronvolts(problem.physical.E_ref)
    scaled_lower=(lower-reference)/scale
    scaled_upper=(upper-reference)/scale
    scaled_step=(scaled_upper-scaled_lower)/(nodes-1)
    _state_require(isfinite(scaled_step) && scaled_step>0,path,"invalid scaled energy step")
    # A bound left unchanged in scaled coordinates can round inward on conversion.
    _state_require(next_lower<=max(lower,scaled_lower*scale+reference) &&
        next_upper>=min(upper,scaled_upper*scale+reference),path,"energy window shrinks")
    physical_nodes=_recipe_energy_nodes((next_upper-next_lower)/((upper-lower)/(nodes-1)),path)
    low_preimage=_scaled_energy_preimage(next_lower,scale,reference)
    high_preimage=_scaled_energy_preimage(next_upper,scale,reference)
    ratio_lower=prevfloat(prevfloat(high_preimage[1]-low_preimage[2])/scaled_step)
    ratio_upper=nextfloat(nextfloat(high_preimage[2]-low_preimage[1])/scaled_step)
    minimum_nodes=ratio_lower<=0 ? 3 : _recipe_energy_nodes(ratio_lower,path)
    maximum_nodes=_recipe_energy_nodes(ratio_upper,path)
    _state_require(next_nodes==physical_nodes || minimum_nodes<=next_nodes<=maximum_nodes,
        path,"energy-node count differs from preserved-step ceil and unit roundtrip")
end
function _validate_recorded_coverage(evidence,problem,path)
    lower,upper=evidence.energy_min_scaled,evidence.energy_max_scaled
    required=(evidence.required_min_scaled,evidence.required_max_scaled)
    _state_require(all(isfinite,(lower,upper,required...)) && lower<upper && required[1]<required[2],
        path,"invalid recorded scaled coverage bounds")
    scale,reference=problem.scales.E₀_eV,_electronvolts(problem.physical.E_ref)
    expected=((evidence.E_min_eV-reference)/scale,(evidence.E_max_eV-reference)/scale)
    _state_require(isequal((lower,upper),expected),path,"scaled bounds differ from recorded grid construction")
    if evidence isa DomainMeasuredEvidence
        _state_require(evidence.scba_quality in _SCBA_RESULT_QUALITIES,
            path,"unknown recorded SCBA quality")
    end
    return lower<=required[1] && upper>=required[2]
end
_recorded_energy_bad(evidence::DomainMeasuredEvidence,policy)=
    maximum((evidence.tail_low,evidence.tail_high,evidence.edge_gamma,evidence.edge_spectral))>
        policy.tail_threshold
function _recorded_domain_finishes(evidence::DomainMeasuredEvidence,revision,policy,problem,path)
    covered=_validate_recorded_coverage(evidence,problem,path)
    return (!_recorded_energy_bad(evidence,policy) && covered) ||
        revision>=policy.maximum_expansions || evidence.scba_quality===:invalid
end
function _recorded_domain_target(row,policy,problem,path)
    evidence=row.evidence
    covered=_validate_recorded_coverage(evidence,problem,path)
    _state_require(row.revision<policy.maximum_expansions,path,"expansion budget already exhausted")
    if evidence isa DomainPreflightEvidence
        _state_require(!covered,path,"already-covered preflight cannot schedule a rebuild")
    else
        _state_require(!_recorded_domain_finishes(evidence,row.revision,policy,problem,path),
            path,"owning measured-domain decision finishes without a rebuild")
    end
    # Replay the owner's scaled ceil BEFORE converting bounds back to eV.
    step=(evidence.energy_max_scaled-evidence.energy_min_scaled)/(evidence.NE-1)
    lower=min(evidence.energy_min_scaled,evidence.required_min_scaled)
    upper=max(evidence.energy_max_scaled,evidence.required_max_scaled)
    nodes=_recipe_energy_nodes((upper-lower)/step,path)
    _state_require(nodes<=policy.maximum_energy_nodes,path,"scaled proposal exceeds declared node cap")
    scale,reference=problem.scales.E₀_eV,_electronvolts(problem.physical.E_ref)
    next_lower=(lower*scale+reference)*u"eV"
    next_upper=(upper*scale+reference)*u"eV"
    if evidence isa DomainMeasuredEvidence && _recorded_energy_bad(evidence,policy)
        # The second ceil is the actual measured-growth expression, in physical units.
        old_lower,old_upper=evidence.E_min_eV*u"eV",evidence.E_max_eV*u"eV"
        width=old_upper-old_lower
        next_lower=min(next_lower,old_lower-policy.growth_fraction*width/2)
        next_upper=max(next_upper,old_upper+policy.growth_fraction*width/2)
        old_step=width/(evidence.NE-1)
        nodes=_recipe_energy_nodes(ustrip(Unitful.NoUnits,(next_upper-next_lower)/old_step),path)
        _state_require(nodes<=policy.maximum_energy_nodes,path,"measured proposal exceeds declared node cap")
    end
    return (_electronvolts(next_lower),_electronvolts(next_upper),nodes)
end
function _validate_revision_energy_transition(row,new,policy,problem,path)
    expected=_recorded_domain_target(row,policy,problem,path)
    _state_require(isequal(new,expected),path,"grid differs from recorded owning adaptation proposal")
end

function _validate_problem_recipe(recipe,problem,production,kernels;current=false)
    for name in (:physical,:scales,:scattering,:models)
        _state_require(_exact_digest(getfield(recipe,name))==_exact_digest(getfield(problem,name)),
            "state.adaptation.recipe.$name","differs from current physical task")
    end
    _state_require(_exact_digest(recipe.algorithms)==_exact_digest(production.algorithms),
        "state.adaptation.recipe.algorithms","differs")
    _state_require(_exact_digest(recipe.kernel_options)==_exact_digest(kernels),
        "state.adaptation.recipe.kernel_options","differs")
    if current
        _state_require(_exact_digest(recipe.numerical)==_exact_digest(problem.numerical),
            "state.adaptation.current_problem_recipe.numerical","differs")
    else
        for name in fieldnames(NumericalParameters)
            name in (:E_min,:E_max,:N_E) && continue
            _state_require(_exact_digest(getfield(recipe.numerical,name))==
                _exact_digest(getfield(problem.numerical,name)),
                "state.adaptation.next_problem_recipe.numerical.$name","undeclared adaptation change")
        end
        old=problem.numerical
        new=recipe.numerical
        _validate_recipe_energy_transition((_electronvolts(old.E_min),_electronvolts(old.E_max),old.N_E),
            (_electronvolts(new.E_min),_electronvolts(new.E_max),new.N_E),problem,
            "state.adaptation.next_problem_recipe.numerical")
    end
end

function validate_exact_state(problem::NEGFProblem,state::ExactSolverState;
    options::SolverOptions,production_options::ProductionOptions,
    execution::ExactExecutionIdentity,domain_adaptation::DomainAdaptationPolicy,
    kernel_options::Union{Nothing,ProductionKernelOptions})
    # Identity rejection precedes array traversal, caches, numerical work and callbacks.
    expected=exact_restart_identity(problem,options,production_options,execution;
        domain_adaptation,kernel_options)
    _state_require(state.identity==expected,"state.identity","problem/policy/execution mismatch")
    _state_require(state.boundary in (:scba_accepted,:outer_accepted,:domain_final,:final),
        "state.boundary","unknown boundary")
    _state_require(state.integrity==_exact_digest(_state_payload(state)),"state.integrity","payload changed")
    _state_require(state.measurement_context.algorithms==production_options.algorithms,
        "state.measurement_context.algorithms","differs from bound policy")
    _state_require(state.measurement_context.backend in (:production,:reference,:not_executed),
        "state.measurement_context.backend","unknown backend")
    _state_require(state.measurement_context.final_state_alignment in (:producer_aligned,:restored_aligned),
        "state.measurement_context.final_state_alignment","unknown alignment provenance")
    n=problem.numerical
    shape=(n.N_E,n.N_k,n.N_b,n.N_b)
    finite=state.boundary ∉ (:domain_final,:final) ||
        (state.final_evidence!==nothing && state.final_evidence.original_converged)
    _state_array(state.published_U,(n.N_z,),"state.published_U";finite)
    _state_array(state.published_n,(n.N_z,),"state.published_n";finite)
    for field in (:Gᴿ,:Gˡ,:Gᵍ,:A)
        _state_array(getfield(state.scba.green,field),shape,"state.scba.green.$field";finite)
    end
    for field in (:condition_number,:dyson_scale)
        _state_array(getfield(state.scba.green,field),(n.N_E,n.N_k),"state.scba.green.$field";finite)
    end
    _state_require(Set(keys(state.scba.scattering))==Set(problem.kernels.enabled),
        "state.scba.scattering","missing or extra mechanisms")
    order=Tuple((problem.kernels.enabled...,:embedding_plus,:embedding_minus))
    _state_require(state.family_order==order,"state.family_order","mixer family order differs")
    for (name,family) in state.scba.scattering
        _state_family(family,shape,"state.scba.scattering.$name";finite)
    end
    for field in (:embedding,:embedding_plus,:embedding_minus)
        _state_family(getfield(state.scba,field),shape,"state.scba.$field";finite)
    end
    mixer=state.scba.mixer_state
    _state_require(mixer.method==production_options.algorithms.mixing,"state.mixer.method","differs")
    _state_require(length(mixer.states)==length(mixer.residuals),"state.mixer.history","length differs")
    _state_require(mixer.method!==:linear || isempty(mixer.states),
        "state.mixer.history","linear mixer cannot carry Anderson history")
    _state_require(length(mixer.states)<=production_options.algorithms.anderson_history_depth,
        "state.mixer.history","exceeds declared depth")
    for (kind,rows) in ((:states,mixer.states),(:residuals,mixer.residuals))
        for (i,row) in enumerate(rows)
            _state_require(length(row)==length(order),"state.mixer.$kind[$i]","missing family")
            for (j,family) in enumerate(row)
                _state_family(family,shape,"state.mixer.$kind[$i][$j]";finite)
            end
        end
    end
    _state_require(state.scba.restart_contract!==nothing,"state.scba.restart_contract","missing")
    _check_exact_restart_contract(state.scba.restart_contract,options,production_options.algorithms)
    inner=state.scba.history
    outer=state.outer_history
    _state_require(all(i->inner[i].ν==i,eachindex(inner)),"state.scba.history","non-contiguous counters")
    _state_require(all(i->outer[i].μ==i,eachindex(outer)),"state.outer_history","non-contiguous counters")
    for (i,row) in enumerate(outer)
        _state_array(row.populations,(n.N_b,),"state.outer_history[$i].populations";finite)
        _state_array(row.density,(n.N_z,),"state.outer_history[$i].density";finite)
    end
    p=state.progress
    _state_require(all(name->getfield(p,name)>=0,fieldnames(SolverProgress)),"state.progress","negative counter")
    _state_require(p.accepted_scba==length(inner)<=options.max_scba,
        "state.progress.accepted_scba","history/cap mismatch")
    _validate_mixer_position(state.mixer_position,mixer,p.accepted_scba,production_options.algorithms)
    _state_require(p.completed_outer==length(outer)<=options.max_poisson,
        "state.progress.completed_outer","history/cap mismatch")
    _state_require(p.outer_strict_streak<=length(outer),"state.progress.outer_strict_streak","exceeds history")
    if p.outer_strict_streak>0
        _state_require(p.outer_strict_streak<=_trailing_pass_count(outer,
            row->poisson_convergence_assessment(row,options.tolerances,options.convergence).passed),
            "state.progress.outer_strict_streak","inconsistent with raw outer history")
    end
    _state_require(1<=p.current_outer<=length(outer)+1,
        "state.progress.current_outer","no owning inner stage")
    if state.boundary in (:domain_final,:final)
        _state_require(p.current_outer in (length(outer),length(outer)+1),
            "state.progress.current_outer","terminal inner stage is not the latest stage")
        _state_require((p.current_outer==length(outer))==scba_accepted(state.scba),
            "state.scba.classification","terminal stage/outer history disagree on accepted inner result")
    end
    before=outer[1:min(p.current_outer-1,length(outer))]
    effective=inner_working_policy(options,before).options
    streak,diagnostic=_history_streaks(state.scba,effective)
    _state_require(p.inner_strict_streak==streak,
        "state.progress.inner_strict_streak","differs from owning inner history")
    _state_require(p.diagnostic_streak==diagnostic,
        "state.progress.diagnostic_streak","differs from owning inner history")
    decision=scba_iteration_decision(inner,effective,streak)
    if state.boundary!==:scba_accepted
        if state.scba.status in (:invalid_candidate,:nonfinite_metrics)
            _state_require(!state.scba.converged && state.scba.quality===:invalid && !isempty(inner) &&
                any(name->!isfinite(getfield(last(inner),name)),
                    (:r_D,:r_A,:r_K,:r_Σ,:r_λ,:r_PSD,:r_caus,:r_roundoff)),
                "state.scba.classification","failure contradicts measured inner history")
        else
            _state_require(decision.terminal && all(name->getfield(state.scba,name)==getfield(decision,name),
                (:converged,:status,:quality)),
                "state.scba.classification","differs from owning inner terminal decision")
        end
        if p.current_outer==length(outer) && p.outer_strict_streak>0
            _state_require(state.scba.converged,"state.progress.outer_strict_streak",
                "current completed outer stage has no strict inner result")
        end
    end
    if state.boundary===:scba_accepted
        c=state.continuation
        _state_require(c isa SCBAContinuation,"state.continuation","SCBA continuation missing")
        _state_require(p.current_outer==length(outer)+1 && !isempty(inner),"state.progress","wrong SCBA boundary")
        _state_require(!state.scba.converged && state.scba.status===:running,
            "state.scba.classification","accepted pre-return snapshot is not running")
        _state_require(all(name->getfield(c.decision,name)==getfield(decision,name),
            (:terminal,:converged,:status,:quality)),"state.continuation.decision","differs from policy/history")
        _state_array(c.next_previous_populations,(n.N_b,),"state.continuation.next_previous_populations")
        _state_require(isequal(c.next_previous_populations,
            _production_state_populations(problem,state.scba.green.Gˡ)),
            "state.continuation.next_previous_populations","differs from accepted Green accumulation")
        _state_require(isfinite(c.next_previous_current),"state.continuation.next_previous_current","non-finite")
        _state_require(isequal(c.next_previous_current,last(inner).J),
            "state.continuation.next_previous_current","differs from accepted row")
        retention=_SCBAWitnessRetention(inner)
        _state_require(c.witness_retention==(retention.first,retention.worst,retention.last),
            "state.continuation.witness_retention","differs from retained history")
        _state_require(c.decision.terminal || c.candidate!==nothing,"state.continuation.candidate","missing")
        if c.candidate!==nothing
            candidate=c.candidate
            _state_require(Set(keys(candidate.scattering))==Set(problem.kernels.enabled),
                "state.continuation.candidate.scattering","missing or extra mechanisms")
            for family in (values(candidate.scattering)...,candidate.embedding,candidate.plus,candidate.minus)
                _state_family(family,shape,"state.continuation.candidate")
            end
        end
    elseif state.boundary===:outer_accepted
        c=state.continuation
        _state_require(c isa OuterContinuation,"state.continuation","outer continuation missing")
        _state_require(p.current_outer==length(outer)>0,"state.progress","wrong outer boundary")
        _state_array(c.next_U,(n.N_z,),"state.continuation.next_U")
        done=p.outer_strict_streak>=options.convergence.required_consecutive_poisson_passes ||
            length(outer)>=options.max_poisson
        _state_require(c.next_action===(done ? :finalize : :next_outer),
            "state.continuation.next_action","differs from policy/counters")
        _state_require(scba_accepted(state.scba),"state.scba","outer boundary has no accepted inner result")
    else
        _state_require(state.continuation===nothing,"state.continuation","terminal state has pending numerical step")
        _state_require(state.final_evidence!==nothing,"state.final_evidence","missing")
        _state_require(state.final_evidence.original_status==state.status,
            "state.final_evidence.original_status","differs from state status")
        if p.current_outer==length(outer)+1
            _state_require(state.status==Symbol(:scba_,state.scba.status),
                "state.status","differs from the unfinished outer stage's inner termination")
        end
        # Reassess the original metrics with the owning band gates, without evaluating
        # a fresh raw map. A passed approximate report is not strict convergence.
        evidence=state.final_evidence
        provisional=NEGFSolution(problem,options,state.published_U,state.published_n,state.scba,
            outer,Dict{Symbol,Any}(),evidence.report,false,state.status)
        strict_report=_validate_solution_metrics(provisional,evidence.report.metrics)
        outer_converged=p.outer_strict_streak>=options.convergence.required_consecutive_poisson_passes
        historically_converged=outer_converged && state.scba.converged && strict_report.passed
        historically_approximate=false
        if !historically_converged && outer_converged && scba_accepted(state.scba) &&
            options.convergence.mode===:adaptive_working && options.convergence.diagnostic_quality.enabled
            approximate=NEGFSolution(problem,_approximate_options(options),state.published_U,
                state.published_n,state.scba,outer,provisional.observables,evidence.report,false,state.status)
            historically_approximate=_validate_solution_metrics(approximate,evidence.report.metrics;
                accept_approximate_scba=true).passed
        end
        _state_require(evidence.report.passed==(strict_report.passed || historically_approximate),
            "state.final_evidence.report.passed","differs from original metrics under the declared quality band")
        _state_require(state.final_evidence.original_converged==historically_converged &&
            (state.status===:converged)==historically_converged &&
            (state.status===:approximate)==historically_approximate,
            "state.final_evidence.original_converged","differs from original history/report classification")
    end
    if domain_adaptation.mode===:none
        _state_require(state.adaptation===nothing,"state.adaptation","unexpected active adaptation")
        _state_require(state.boundary!==:domain_final,"state.boundary","domain final without adaptive wrapper")
    else
        a=state.adaptation
        _state_require(a!==nothing,"state.adaptation","missing")
        _state_require(a.policy==domain_adaptation && a.kernel_options==kernel_options,
            "state.adaptation.policy","rebuild policy differs")
        _state_require(a.current_problem==state.identity.problem,"state.adaptation.current_problem","differs")
        _validate_problem_recipe(a.current_problem_recipe,problem,production_options,kernel_options;current=true)
        _state_require(0<=a.expansions<=domain_adaptation.maximum_expansions,"state.adaptation.expansions","outside declared cap")
        _state_require(a.expansions==0 || n.N_E<=domain_adaptation.maximum_energy_nodes,
            "state.adaptation.current_problem_recipe.numerical.N_E","expanded domain exceeds declared cap")
        _state_require(isfinite(a.elapsed_domain_seconds) && a.elapsed_domain_seconds>=0,
            "state.adaptation.elapsed_domain_seconds","missing or invalid recorded domain cost")
        _state_require(state.boundary!==:final || a.next_action===:finished,
            "state.adaptation.next_action","global final has pending action")
        _state_require(state.boundary!==:domain_final || a.next_action!==:finished,
            "state.adaptation.next_action","domain final is already global final")
        _state_require(a.next_action in (:assess_measured_domain,:cold_rebuild,:finished),
            "state.adaptation.next_action","unknown checkpoint action")
        _state_require(a.next_action!==:cold_rebuild || state.boundary===:domain_final,
            "state.adaptation.next_action","cold rebuild requires a completed domain")
        _state_require(a.next_action!==:finished || state.boundary===:final,
            "state.adaptation.next_action","finished requires global final")
        committed=a.next_action in (:cold_rebuild,:finished)
        _state_require(a.decision_committed==committed,
            "state.adaptation.decision_committed","differs from action")
        _state_require(a.next_action===:cold_rebuild || a.next_problem_recipe===nothing,
            "state.adaptation.next_problem_recipe","unexpected pending recipe")
        _state_require(a.decision_id==length(a.revisions),"state.adaptation.decision_id","differs from committed revision count")
        indices=[row.revision for row in a.revisions]
        expected_indices=collect(0:(a.expansions-(committed ? 0 : 1)))
        _state_require(indices==expected_indices,"state.adaptation.revisions","missing revision or invalid order/counters")
        previous_revision=nothing
        for row in a.revisions
            _state_require((row.source===:preflight_estimate && row.evidence isa DomainPreflightEvidence) ||
                (row.source===:measured_state && row.evidence isa DomainMeasuredEvidence),
                "state.adaptation.revisions[$(row.revision)].source","differs from typed evidence")
            grid=(row.evidence.E_min_eV,row.evidence.E_max_eV,row.evidence.NE)
            path="state.adaptation.revisions[$(row.revision)].grid"
            _recipe_energy_window(grid,problem,path)
            _validate_recorded_coverage(row.evidence,problem,path)
            _state_require(row.revision==0 || row.evidence.NE<=domain_adaptation.maximum_energy_nodes,
                path,"expanded domain exceeds declared cap")
            previous_revision===nothing || _validate_revision_energy_transition(previous_revision,grid,
                domain_adaptation,problem,path)
            previous_revision=row
            if row.evidence isa DomainMeasuredEvidence
                _state_require(isfinite(row.evidence.seconds) && row.evidence.seconds>=0,
                    "state.adaptation.revisions[$(row.revision)].seconds","invalid measured cost")
                row.revision==a.expansions && _state_require(row.evidence.seconds<=a.elapsed_domain_seconds,
                    "state.adaptation.elapsed_domain_seconds","precedes committed measured receipt")
            end
            if row.revision==a.expansions
                _state_require(row.evidence isa DomainMeasuredEvidence,
                    "state.adaptation.revisions[$(row.revision)]","completed domain has no measured receipt")
                evidence=row.evidence
                _state_require(isequal(grid,(_electronvolts(n.E_min),_electronvolts(n.E_max),n.N_E)),
                    path,"differs from current numerical domain")
                g,b=problem.grids,problem.basis
                _state_require(evidence.discretization_id==_representation_digest(
                    g.x,g.wˣ,g.ε,g.wᴱ,g.κ,g.wᵏ,g.qᶻ,g.φ) &&
                    evidence.basis_id==_representation_digest(b.Φ,b.H₀,b.M⁻¹) &&
                    evidence.hartree_id==_representation_digest(state.published_U),
                    "state.adaptation.revisions[$(row.revision)].identity","differs from current arrays")
                _state_require(state.final_evidence!==nothing && evidence.status==String(state.status) &&
                    evidence.scba_quality===state.scba.quality &&
                    evidence.strict_accepted==state.final_evidence.original_converged &&
                    isequal(evidence.metrics,state.final_evidence.report.metrics),
                    "state.adaptation.revisions[$(row.revision)]","contradicts original final evidence")
            end
        end
        if previous_revision!==nothing && previous_revision.revision<a.expansions
            _validate_revision_energy_transition(previous_revision,
                (_electronvolts(n.E_min),_electronvolts(n.E_max),n.N_E),domain_adaptation,problem,
                "state.adaptation.current_problem_recipe.numerical")
        end
        if committed
            finishes=_recorded_domain_finishes(last(a.revisions).evidence,a.expansions,
                domain_adaptation,problem,"state.adaptation.current_measured_receipt")
            _state_require((a.next_action===:finished)==finishes,"state.adaptation.next_action",
                "differs from owning recorded tail/coverage/quality/budget decision")
        end
        if a.next_action===:cold_rebuild
            _state_require(state.scba.quality!==:invalid,"state.adaptation.next_action",
                "invalid SCBA terminates adaptive execution without a cold rebuild")
            _state_require(a.expansions<domain_adaptation.maximum_expansions,
                "state.adaptation.expansions","cold rebuild exceeds declared cap")
            _state_require(a.next_problem_recipe!==nothing,"state.adaptation.next_problem_recipe","missing")
            _validate_problem_recipe(a.next_problem_recipe,problem,production_options,kernel_options)
            next_n=a.next_problem_recipe.numerical
            _validate_revision_energy_transition(last(a.revisions),
                (_electronvolts(next_n.E_min),_electronvolts(next_n.E_max),next_n.N_E),
                domain_adaptation,problem,"state.adaptation.next_problem_recipe.numerical")
            _state_require(a.next_problem_recipe.numerical.N_E<=domain_adaptation.maximum_energy_nodes,
                "state.adaptation.next_problem_recipe.numerical.N_E","exceeds declared cap")
        end
    end
    return nothing
end

function restore_exact_state(problem::NEGFProblem,state::ExactSolverState;
    options::SolverOptions,production_options::ProductionOptions,
    execution::ExactExecutionIdentity,domain_adaptation::DomainAdaptationPolicy,
    kernel_options::Union{Nothing,ProductionKernelOptions})
    validate_exact_state(problem,state;options,production_options,execution,domain_adaptation,kernel_options)
    return RestoredSolverState(problem,state,options,production_options,execution,
        domain_adaptation,kernel_options,Val(:validated))
end
