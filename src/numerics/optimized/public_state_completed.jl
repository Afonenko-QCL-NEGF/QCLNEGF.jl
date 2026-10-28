struct RevalidationReceipt
    status::Symbol
    reason::Union{Nothing,Symbol}
    detail::String
    context::AcceptanceContext
    report::ConvergenceReport
end
function RevalidationReceipt(;status,reason,detail,context,report)
    status in (:measured,:error) || throw(ArgumentError("revalidation.status is unknown"))
    status===:error && reason===nothing && throw(ArgumentError("revalidation.reason is missing"))
    RevalidationReceipt(status,reason,String(detail),context,deepcopy(report))
end
struct IndependentOracleReceipt
    status::Symbol
    reason::Union{Nothing,Symbol}
    detail::String
    context::Union{Nothing,AcceptanceContext}
    report::Union{Nothing,ConvergenceReport}
end
function IndependentOracleReceipt(;status,reason,detail,context,report)
    status in (:measured,:not_measured,:unsupported,:error) ||
        throw(ArgumentError("oracle.status is unknown"))
    status!==:measured && reason===nothing && throw(ArgumentError("oracle.reason is missing"))
    IndependentOracleReceipt(status,reason,String(detail),context,deepcopy(report))
end
struct CompletedStateEvidence
    original::FinalEvidenceSnapshot
    revalidation::RevalidationReceipt
    independent_oracle::IndependentOracleReceipt
end
CompletedStateEvidence(;original,revalidation,independent_oracle)=
    CompletedStateEvidence(deepcopy(original),deepcopy(revalidation),deepcopy(independent_oracle))

function _terminal_arrays_finite(state)
    all(isfinite,state.published_U) && all(isfinite,state.published_n) || return false
    for field in fieldnames(GreenState)
        all(isfinite,getfield(state.scba.green,field)) || return false
    end
    for family in (values(state.scba.scattering)...,state.scba.embedding,
                   state.scba.embedding_plus,state.scba.embedding_minus)
        all(component->all(isfinite,getfield(family,component)),(:Σᴿ,:Σˡ,:Σᵍ)) || return false
    end
    return true
end
function _unmeasurable_completed(restart,state,reason,detail,oracle_request)
    observations=Dict{Symbol,Any}(Symbol(first(pair))=>_core_plain(last(pair))
        for pair in state.final_evidence.observations.entries if first(pair)!="final_fixed_point_audit")
    observations[:warnings]=Dict{String,Any}[_core_plain(row) for row in state.warnings]
    observations[:acceptance_measurement_context]=Dict{String,Any}(
        "backend"=>"not_executed","hilbert_roundoff_measured"=>false,
        "map_kernels_shared_with_solver"=>nothing,"final_state_alignment"=>"restored_aligned")
    report=ConvergenceReport(false,Dict{Symbol,Float64}(),[detail])
    context=AcceptanceContext(backend=:not_executed,algorithms=restart.production_options.algorithms,
        hilbert_roundoff_measured=nothing,map_kernels_shared_with_solver=nothing,
        final_state_alignment=:restored_aligned)
    selected=RevalidationReceipt(:error,reason,detail,context,report)
    independent=oracle_request===:none ?
        IndependentOracleReceipt(:not_measured,:not_requested,"",nothing,nothing) :
        IndependentOracleReceipt(:error,reason,detail,nothing,nothing)
    observations[:restoration_evidence]=CompletedStateEvidence(deepcopy(state.final_evidence),selected,independent)
    result=NEGFSolution(restart.problem,restart.options,state.published_U,state.published_n,state.scba,
        state.outer_history,observations,report,false,:validation_failed)
    observations[:quality]=solution_quality(result)
    observations[:termination_reason]="validation_failed"
    return result
end

function _completed_candidate(problem,scba,U,production)
    cache=build_production_cache(problem;options=production)
    scattering,roundoff=_scattering_candidate_production(problem,scba.green,cache,production;
        return_roundoff=true)
    embedding,plus,minus=_embedding_self_energy_selected(problem,scba.green,cache,production;
        scattering,h=project_hamiltonians(problem,U))
    return FixedPointAuditCandidate(scattering,embedding,plus,minus,roundoff)
end
function _independent_completed_oracle(provisional,request)
    request===:none && return IndependentOracleReceipt(:not_measured,:not_requested,"",nothing,nothing)
    algorithms=_algorithm_options_from_solution(provisional)
    separate=algorithms.energy_shift===:sparse_plan &&
        provisional.problem.energy_shift_discretization===:nodal_linear
    supported=separate && algorithms.solver_backend===:production && isempty(provisional.problem.kernels.enabled) &&
        algorithms.embedding===:full_cell_resolvent && algorithms.self_energy_structure===:full &&
        algorithms.transverse_momentum===:resolved && algorithms.retarded_real_part===:kramers_kronig
    if !supported
        reason=algorithms.embedding===:finite_chain ? :finite_chain_assembly_unavailable :
            algorithms.energy_shift in (:dense,:conservative_pair) ? :shared_selected_reference_implementation :
            :independent_map_not_qualified
        return IndependentOracleReceipt(:unsupported,reason,
            "No independently qualified equivalent assembly for this saved selected map.",nothing,nothing)
    end
    context=AcceptanceContext(;backend=:reference,algorithms,
        hilbert_roundoff_measured=true,map_kernels_shared_with_solver=false,
        final_state_alignment=:restored_aligned)
    try
        candidate=_reference_audit_candidate(provisional.problem,provisional.scba)
        independent=deepcopy(provisional)
        report=validate_solution(independent;candidate)
        return IndependentOracleReceipt(:measured,nothing,"literal bulk embedding; no scattering",context,report)
    catch error
        error isa InterruptException && rethrow()
        return IndependentOracleReceipt(:error,:independent_measurement_failed,
            sprint(showerror,error),context,nothing)
    end
end

"""Revalidate a global final state on its saved map, without SCBA/Poisson/polish iterations."""
function restore_completed_solution(restart::RestoredSolverState;
    revalidation::Symbol,independent_oracle::Symbol)
    revalidation===:selected_map || throw(ArgumentError("revalidation must be :selected_map"))
    independent_oracle in (:none,:if_supported) || throw(ArgumentError("independent_oracle is unknown"))
    state=_restore_binding(restart.problem,restart)
    _state_require(state.boundary===:final,"restart.boundary","completed restore requires global final")
    problem,options,production=restart.problem,restart.options,restart.production_options
    _terminal_arrays_finite(state) || return _unmeasurable_completed(restart,state,
        :nonfinite_terminal_state,"Stored terminal arrays are non-finite; no numerical revalidation was performed.",
        independent_oracle)
    observations=try
        _collect_observables(problem,state.scba,state.published_U)
    catch error
        error isa InterruptException && rethrow()
        return _unmeasurable_completed(restart,state,:observable_measurement_failed,
            sprint(showerror,error),independent_oracle)
    end
    for pair in state.final_evidence.observations.entries
        observations[Symbol(first(pair))]=_core_plain(last(pair))
    end
    pop!(observations,:final_fixed_point_audit,nothing)
    observations[:warnings]=Dict{String,Any}[_core_plain(row) for row in state.warnings]
    observations[:restart_contract]=_solver_restart_contract(options,production.algorithms)
    roundoff_measured=production.verify_fft_roundoff || all(==(:electron_electron),problem.kernels.enabled)
    context=AcceptanceContext(backend=:production,algorithms=production.algorithms,
        hilbert_roundoff_measured=roundoff_measured,map_kernels_shared_with_solver=true,
        final_state_alignment=:restored_aligned)
    observations[:acceptance_measurement_context]=Dict{String,Any}(
        "backend"=>"production","map_kernels_shared_with_solver"=>true,
        "final_state_alignment"=>"restored_aligned","hilbert_roundoff_measured"=>roundoff_measured,
        "fresh_candidate_required"=>true)
    provisional=NEGFSolution(problem,options,state.published_U,state.published_n,state.scba,
        state.outer_history,observations,deepcopy(state.final_evidence.report),false,state.status)
    candidate=nothing; failure=nothing
    try
        candidate=_completed_candidate(problem,state.scba,state.published_U,production)
    catch error
        error isa InterruptException && rethrow()
        failure=sprint(showerror,error)
    end
    if failure!==nothing
        context=AcceptanceContext(backend=:production,algorithms=production.algorithms,
            hilbert_roundoff_measured=nothing,map_kernels_shared_with_solver=true,
            final_state_alignment=:restored_aligned)
        observations[:acceptance_measurement_context]["hilbert_roundoff_measured"]=false
    end
    outer_converged=state.progress.outer_strict_streak>=options.convergence.required_consecutive_poisson_passes
    quality=try
        _final_quality_report(provisional,outer_converged;candidate,audit_failure=failure)
    catch error
        error isa InterruptException && rethrow()
        failure=sprint(showerror,error)
        (;report=ConvergenceReport(false,Dict{Symbol,Float64}(),
            ["Completed state could not be revalidated: $failure"]),converged=false,approximate=false)
    end
    _publish_charge_diagnostics!(observations,quality.report)
    selected=RevalidationReceipt(failure===nothing ? :measured : :error,
        failure===nothing ? nothing : :selected_map_measurement_failed,
        failure===nothing ? "" : failure,context,deepcopy(quality.report))
    independent=_independent_completed_oracle(provisional,independent_oracle)
    observations[:restoration_evidence]=CompletedStateEvidence(deepcopy(state.final_evidence),selected,independent)
    status=quality.converged ? :converged : quality.approximate ? :approximate :
        quality.report.passed ? state.status : :validation_failed
    result=NEGFSolution(problem,options,state.published_U,state.published_n,state.scba,state.outer_history,
        observations,quality.report,quality.converged,status)
    observations[:quality]=solution_quality(result)
    observations[:termination_reason]=String(status)
    return result
end
