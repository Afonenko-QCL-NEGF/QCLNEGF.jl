module PublicCompletedRestoreTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")

function completed_map_fixture(f,production)
    old=public_state_failed_scba(f.state.scba)
    scba=SCBAResult(old.green,old.scattering,old.embedding,old.embedding_plus,old.embedding_minus,
        old.history,old.converged,old.status,old.quality,
        QCLNEGF._solver_restart_contract(f.options,production.algorithms),old.mixer_state)
    execution=public_state_execution(production)
    identity=exact_restart_identity(f.problem,f.options,production,execution;
        domain_adaptation=f.policy,kernel_options=nothing)
    original=FinalEvidenceSnapshot(report=ConvergenceReport(false,Dict(:r_K=>NaN),["synthetic map fixture"]),
        original_converged=false,original_status=:scba_invalid_candidate,
        observations=CoreRecord(Dict("quality"=>"unconverged")))
    state=ExactSolverState(identity=identity,boundary=:final,published_U=f.state.published_U,
        published_n=f.state.published_n,scba=scba,mixer_position=f.state.mixer_position,family_order=f.state.family_order,
        outer_history=f.state.outer_history,progress=f.state.progress,continuation=nothing,
        status=:scba_invalid_candidate,warnings=f.state.warnings,policy_evaluations=(),adaptation=nothing,
        measurement_context=AcceptanceContext(backend=:production,algorithms=production.algorithms,
            hilbert_roundoff_measured=false,map_kernels_shared_with_solver=true,
            final_state_alignment=:producer_aligned),final_evidence=original)
    handle=restore_exact_state(f.problem,state;options=f.options,production_options=production,
        execution,domain_adaptation=f.policy,kernel_options=nothing)
    return state,handle
end

@testset "Completed restore preserves postmortem arrays and original missing evidence" begin
    @test isdefined(QCLNEGF,:restore_completed_solution)
    if isdefined(QCLNEGF,:restore_completed_solution)
        f=public_state_fixture()
        original=FinalEvidenceSnapshot(report=ConvergenceReport(false,
            Dict(:r_K=>NaN,:r_U=>Inf),["synthetic postmortem"]),
            original_converged=false,original_status=:scba_invalid_candidate,
            observations=CoreRecord(Dict("acceptance_measurement_context"=>Dict(
                "hilbert_roundoff_measured"=>false))))
        state=ExactSolverState(identity=f.state.identity,boundary=:final,
            published_U=f.state.published_U,published_n=f.state.published_n,
            scba=public_state_failed_scba(f.state.scba),mixer_position=f.state.mixer_position,
            family_order=f.state.family_order,outer_history=f.state.outer_history,
            progress=f.state.progress,continuation=nothing,status=:scba_invalid_candidate,
            warnings=f.state.warnings,policy_evaluations=(),adaptation=nothing,
            measurement_context=f.state.measurement_context,final_evidence=original)
        handle=restore_exact_state(f.problem,state;options=f.options,
            production_options=f.production,execution=f.execution,
            domain_adaptation=f.policy,kernel_options=nothing)
        before=deepcopy(state)
        result=restore_completed_solution(handle;revalidation=:selected_map,
            independent_oracle=:if_supported)
        @test result isa NEGFSolution
        @test result.Uᴴ == before.published_U
        @test result.n == before.published_n
        @test result.scba.green.Gᴿ == before.scba.green.Gᴿ
        @test state.integrity == before.integrity
        receipt=result.observables[:restoration_evidence]
        @test isnan(receipt.original.report.metrics[:r_K])
        @test isinf(receipt.original.report.metrics[:r_U])
        @test receipt.original.observations["acceptance_measurement_context"]["hilbert_roundoff_measured"].value === false
        @test receipt.revalidation.context.algorithms == f.production.algorithms
        @test receipt.independent_oracle.status in (:measured,:unsupported,:error)
        @test receipt.revalidation.status===:measured
        @test receipt.revalidation.reason===nothing
        @test !receipt.revalidation.report.passed
        @test !solution_scientific_assessment(result)["scientific_accepted"]
        result.scba.green.Gᴿ[1]=7.0+8im
        @test state.scba.green.Gᴿ[1] == before.scba.green.Gᴿ[1]
        @test_throws ArgumentError solve_production(f.problem;restart=handle)
    end
end

@testset "A shared dense embedding implementation is not labelled an independent oracle" begin
    f=public_state_fixture()
    production=QCLNEGF.with_production_options(f.production;
        algorithms=AlgorithmOptions(energy_shift=:dense))
    state,handle=completed_map_fixture(f,production)
    result=restore_completed_solution(handle;revalidation=:selected_map,independent_oracle=:if_supported)
    receipt=result.observables[:restoration_evidence]
    @test receipt.independent_oracle.status===:unsupported
    @test receipt.independent_oracle.reason===:shared_selected_reference_implementation
end

@testset "Finite-chain revalidation does not silently substitute the legacy bulk closure" begin
    f=public_state_fixture()
    production=QCLNEGF.with_production_options(f.production;
        algorithms=AlgorithmOptions(embedding=:finite_chain))
    state,handle=completed_map_fixture(f,production)
    N=QCLNEGF.QCLNumerics
    selected=N._completed_candidate(f.problem,state.scba,state.published_U,production)
    bulk=N._reference_audit_candidate(f.problem,state.scba)
    @test norm(selected.embedding.Σᴿ-bulk.embedding.Σᴿ)>0
    result=restore_completed_solution(handle;revalidation=:selected_map,independent_oracle=:if_supported)
    receipt=result.observables[:restoration_evidence]
    @test receipt.revalidation.context.algorithms.embedding===:finite_chain
    @test receipt.independent_oracle.status===:unsupported
    @test receipt.independent_oracle.reason===:finite_chain_assembly_unavailable
    @test !solution_scientific_assessment(result)["scientific_accepted"]
end

@testset "An undefined attempted map does not claim an unexecuted roundoff measurement" begin
    # Manufactured invalid postmortem Green pair, not an accepted physical state or a solve.
    f=public_state_fixture(scattering_options=ScatteringOptions(
        LO=true,acoustic=false,impurity=false,IFR=false,alloy=false),
        physical_models=PhysicalModelOptions(lo_population=:rate_balance,
            lo_decay_ps=1.0,lo_mode_density_per_m3=1e24))
    f.state.scba.green.Gᵍ .= 1.5im
    _,handle=completed_map_fixture(f,f.production)
    result=restore_completed_solution(handle;revalidation=:selected_map,independent_oracle=:none)
    receipt=result.observables[:restoration_evidence]
    @test receipt.revalidation.status===:error
    @test occursin("negative electron-phonon transition rate",receipt.revalidation.detail)
    @test receipt.revalidation.context.hilbert_roundoff_measured===nothing
    @test !solution_scientific_assessment(result)["scientific_accepted"]
end

@testset "Non-finite terminal data remain readable without a fabricated fresh audit" begin
    f=public_state_fixture()
    source,_=completed_map_fixture(f,f.production)
    U=fill(0.125,length(source.published_U)); U[1]=NaN
    state=ExactSolverState(identity=source.identity,boundary=:final,published_U=U,
        published_n=source.published_n,scba=source.scba,mixer_position=source.mixer_position,family_order=source.family_order,
        outer_history=source.outer_history,progress=source.progress,continuation=nothing,
        status=source.status,warnings=source.warnings,policy_evaluations=(),adaptation=nothing,
        measurement_context=source.measurement_context,final_evidence=source.final_evidence)
    handle=restore_exact_state(f.problem,state;options=f.options,production_options=f.production,
        execution=f.execution,domain_adaptation=f.policy,kernel_options=nothing)
    result=try
        restore_completed_solution(handle;revalidation=:selected_map,independent_oracle=:none)
    catch error
        error
    end
    @test typeof(result)===NEGFSolution
    if result isa NEGFSolution
        @test isnan(result.Uᴴ[1])
        @test result.Uᴴ[2:end]==fill(0.125,length(U)-1)
        @test result.observables[:restoration_evidence].revalidation.status===:error
        @test result.observables[:restoration_evidence].revalidation.reason===:nonfinite_terminal_state
        @test solution_scientific_assessment(result)["final_fixed_point_audit"]===nothing
        @test !solution_scientific_assessment(result)["scientific_accepted"]
    end
end

@testset "Historical strict and approximate final bands remain distinct on restore" begin
    # The actual final-quality evaluator owns this branch. The altered raw candidate and
    # permissive thresholds control classification only; this is not physics evidence.
    N=QCLNEGF.QCLNumerics
    problem=public_state_problem()
    tolerances=SolverTolerances(; (name=>1e100 for name in fieldnames(SolverTolerances))...)
    strict_options=SolverOptions(α_Σ=0.1,max_scba=4,max_poisson=3,tolerances=tolerances,
        convergence=ConvergencePolicy(minimum_scba_iterations=4,
            required_consecutive_scba_passes=1,required_consecutive_poisson_passes=1,
            stagnation_window=0,stagnation_relative_improvement=0.0))
    production=ProductionOptions(worker_count=1,checkpoint_every_scba=1,checkpoint_every_outer=1,
        progress_every_scba=0,progress_every_outer=0)
    execution=public_state_execution(production)
    states=ExactSolverState[]
    result=solve_production(problem;options=strict_options,production_options=production,
        execution,checkpoint_sink=state->push!(states,state))
    strict_state=only(filter(state->state.boundary===:final,states))
    @test strict_state.final_evidence.original_converged

    quality_policy=DiagnosticQualityPolicy(enabled=true,keldysh_threshold=1e100,
        self_energy_threshold=1e100,normalization_threshold=1e100,observable_threshold=1e100,
        positivity_threshold=1e100,causality_threshold=1e100,outer_threshold=1e100)
    policy=ConvergencePolicy(; (name=>(name===:mode ? :adaptive_working :
        name===:diagnostic_quality ? quality_policy : getfield(strict_options.convergence,name))
        for name in fieldnames(ConvergencePolicy))...)
    band_tolerances=SolverTolerances(; (name=>(name===:r_K ?
        nextfloat(last(strict_state.scba.history).r_K) : getfield(tolerances,name))
        for name in fieldnames(SolverTolerances))...)
    options=SolverOptions(; (name=>(name===:convergence ? policy : name===:tolerances ?
        band_tolerances : getfield(strict_options,name)) for name in fieldnames(SolverOptions))...)

    function rebind(state,options;domain_adaptation=DomainAdaptationPolicy(),kernel_options=nothing,
        changes...)
        old=state.scba
        scba=SCBAResult(old.green,old.scattering,old.embedding,old.embedding_plus,
            old.embedding_minus,old.history,old.converged,old.status,old.quality,
            QCLNEGF._solver_restart_contract(options,production.algorithms),old.mixer_state)
        p=state.progress
        effective=N.inner_working_policy(options,state.outer_history[1:p.current_outer-1]).options
        strict,diagnostic=N._history_streaks(scba,effective)
        identity=exact_restart_identity(problem,options,production,execution;
            domain_adaptation,kernel_options)
        return replace_public_state(state;identity,scba,progress=SolverProgress(p.completed_outer,
            p.current_outer,p.accepted_scba,p.outer_strict_streak,strict,diagnostic),changes...)
    end
    rebound=rebind(strict_state,options)
    provisional=NEGFSolution(problem,options,rebound.published_U,rebound.published_n,rebound.scba,
        rebound.outer_history,deepcopy(result.observables),rebound.final_evidence.report,false,:running)
    candidate=N._completed_candidate(problem,rebound.scba,rebound.published_U,production)
    candidate.embedding.Σˡ .*= 1e9
    quality=N._final_quality_report(provisional,true;candidate)
    @test quality.approximate && quality.report.passed && !quality.converged
    @test rebound.scba.converged
    @test N.scba_iteration_decision(rebound.scba.history,options,
        rebound.progress.inner_strict_streak).converged
    original=FinalEvidenceSnapshot(report=quality.report,original_converged=false,
        original_status=:approximate,observations=CoreRecord(Dict("quality"=>"approximate")))
    approximate=replace_public_state(rebound;status=:approximate,final_evidence=original)
    @test validate_exact_state(problem,approximate;options,production_options=production,
        execution,domain_adaptation=DomainAdaptationPolicy(),kernel_options=nothing)===nothing

    # Neither a false strict flag nor an approximate label can hide strict-passing metrics.
    for status in (:validation_failed,:approximate)
        wrong=replace_public_state(strict_state;status,final_evidence=FinalEvidenceSnapshot(
            report=strict_state.final_evidence.report,original_converged=false,
            original_status=status,observations=original.observations))
        @test_throws ArgumentError validate_exact_state(problem,wrong;options=strict_options,
            production_options=production,execution,domain_adaptation=DomainAdaptationPolicy(),
            kernel_options=nothing)
    end
    false_strict_report=ConvergenceReport(false,strict_state.final_evidence.report.metrics,
        ["corrupted strict passed flag"])
    wrong=replace_public_state(strict_state;status=:validation_failed,
        final_evidence=FinalEvidenceSnapshot(report=false_strict_report,original_converged=false,
            original_status=:validation_failed,observations=original.observations))
    @test_throws ArgumentError validate_exact_state(problem,wrong;options=strict_options,
        production_options=production,execution,domain_adaptation=DomainAdaptationPolicy(),
        kernel_options=nothing)
    wrong=replace_public_state(approximate;status=:converged,final_evidence=FinalEvidenceSnapshot(
        report=quality.report,original_converged=true,original_status=:converged,
        observations=original.observations))
    @test_throws ArgumentError validate_exact_state(problem,wrong;options,
        production_options=production,execution,domain_adaptation=DomainAdaptationPolicy(),
        kernel_options=nothing)
    for mode in (:strict_fail_fast,:adaptive_working)
        q=DiagnosticQualityPolicy(; (name=>(name===:enabled ? false :
            getfield(quality_policy,name)) for name in fieldnames(DiagnosticQualityPolicy))...)
        disabled=ConvergencePolicy(; (name=>(name===:mode ? mode :
            name===:diagnostic_quality ? (mode===:adaptive_working ? q : quality_policy) :
            getfield(policy,name)) for name in fieldnames(ConvergencePolicy))...)
        disabled_options=SolverOptions(; (name=>(name===:convergence ? disabled :
            getfield(options,name)) for name in fieldnames(SolverOptions))...)
        disabled_state=rebind(approximate,disabled_options)
        @test disabled_state.scba.converged
        @test_throws ArgumentError validate_exact_state(problem,disabled_state;
            options=disabled_options,production_options=production,execution,
            domain_adaptation=DomainAdaptationPolicy(),kernel_options=nothing)
    end
    bad_report=ConvergenceReport(false,quality.report.metrics,["corrupted approximate passed flag"])
    wrong=replace_public_state(approximate;final_evidence=FinalEvidenceSnapshot(report=bad_report,
        original_converged=false,original_status=:approximate,observations=original.observations))
    @test_throws ArgumentError validate_exact_state(problem,wrong;options,
        production_options=production,execution,domain_adaptation=DomainAdaptationPolicy(),
        kernel_options=nothing)

    handle=restore_exact_state(problem,approximate;options,production_options=production,
        execution,domain_adaptation=DomainAdaptationPolicy(),kernel_options=nothing)
    fresh=restore_completed_solution(handle;revalidation=:selected_map,independent_oracle=:none)
    receipt=fresh.observables[:restoration_evidence]
    @test !receipt.original.original_converged
    @test receipt.original.original_status===:approximate
    @test receipt.original.report.passed
    @test isequal(receipt.original.report.metrics,quality.report.metrics)
    @test receipt.revalidation.report.metrics[:r_K] != quality.report.metrics[:r_K]
    @test isequal(approximate.final_evidence.report.metrics,quality.report.metrics)

    domain_policy=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=0,
        maximum_energy_nodes=31)
    kernels=ProductionKernelOptions()
    identity=exact_restart_identity(problem,options,production,execution;
        domain_adaptation=domain_policy,kernel_options=kernels)
    adaptation=DomainAdaptationState(policy=domain_policy,kernel_options=kernels,expansions=0,
        elapsed_domain_seconds=0.125,current_problem=identity.problem,
        current_problem_recipe=N._problem_recipe(problem,production.algorithms,kernels),
        revisions=(),next_action=:assess_measured_domain,decision_id=0,decision_committed=false,
        next_problem_recipe=nothing)
    domain=rebind(approximate,options;domain_adaptation=domain_policy,kernel_options=kernels,
        boundary=:domain_final,adaptation)
    @test validate_exact_state(problem,domain;options,production_options=production,execution,
        domain_adaptation=domain_policy,kernel_options=kernels)===nothing
    handle=restore_exact_state(problem,domain;options,production_options=production,execution,
        domain_adaptation=domain_policy,kernel_options=kernels)
    resumed_states=ExactSolverState[]
    rows=ScientificHistorySnapshot[]
    resumed=solve_adaptive_production(problem;restart=handle,
        checkpoint_sink=state->push!(resumed_states,state),history_observer=row->push!(rows,row))
    @test isempty(rows)
    @test !resumed.converged && resumed.status===:approximate && resumed.report.passed
    @test resumed.scba.green.Gᴿ==domain.scba.green.Gᴿ
    @test isequal(resumed.report.metrics,domain.final_evidence.report.metrics)
    global_final=last(resumed_states)
    @test global_final.boundary===:final
    @test !only(global_final.adaptation.revisions).evidence.strict_accepted
    @test validate_exact_state(problem,global_final;options,production_options=production,execution,
        domain_adaptation=domain_policy,kernel_options=kernels)===nothing
end
end
