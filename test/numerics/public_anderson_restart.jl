module PublicAndersonRestartTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")
const N=QCLNEGF.QCLNumerics
function causal_family(gamma)
    SelfEnergyFamily(fill(ComplexF64(-im*gamma),1,1,1,1),
        fill(ComplexF64(im*gamma),1,1,1,1),fill(ComplexF64(-im*gamma),1,1,1,1))
end

@testset "Saved Anderson history matches a literal two-channel constrained least-squares oracle" begin
    algorithms=AlgorithmOptions(mixing=:anderson,anderson_history_depth=2,
        anderson_regularization=0.0,anderson_damping=1.0)
    production=ProductionOptions(worker_count=1,algorithms=algorithms)
    workspace=N._AndersonWorkspace()
    state=[causal_family(2.0),causal_family(4.0)]
    N._anderson_mix!(workspace,state,[causal_family(3.0),causal_family(6.0)],algorithms,1.0,production)
    @test real(-im*state[1].Σˡ[1]) == 3.0
    N._anderson_mix!(workspace,state,[causal_family(5.0),causal_family(7.0)],algorithms,1.0,production)
    # Gram = 3*[5 4;4 5]; constrained weights are exactly (1/2,1/2).
    @test real(-im*state[1].Σˡ[1]) ≈ 4.0 atol=1e-13
    @test real(-im*state[2].Σˡ[1]) ≈ 6.5 atol=1e-13
    saved=SCBAMixerState(:anderson,deepcopy(workspace.states),deepcopy(workspace.residuals))
    resumed=N._AndersonWorkspace(deepcopy(saved.states),deepcopy(saved.residuals))
    continued=deepcopy(state)
    candidate=[causal_family(5.0),causal_family(6.75)]
    N._anderson_mix!(workspace,state,candidate,algorithms,1.0,production)
    N._anderson_mix!(resumed,continued,candidate,algorithms,1.0,production)
    # Oldest row evicted: Gram = [15 27/4;27/4 51/16].
    # Independent 2x2 inverse gives weights (-19/25,44/25), hence (5,164/25).
    @test real(-im*continued[1].Σˡ[1]) ≈ 5.0 atol=1e-12
    @test real(-im*continued[2].Σˡ[1]) ≈ 164/25 atol=1e-12
    @test length(resumed.states)==2
    for j=1:2, component in (:Σᴿ,:Σˡ,:Σᵍ)
        @test getfield(continued[j],component) == getfield(state[j],component)
    end
    @test real(-im*resumed.residuals[1][1].Σˡ[1]) ≈ 2.0 atol=1e-13
    @test real(-im*resumed.residuals[2][2].Σˡ[1]) ≈ 0.25 atol=1e-13
end

@testset "Public continuation restores nonempty Anderson payload and raw history" begin
    problem=public_state_problem()
    tolerances=SolverTolerances(; (field=>1e100 for field in fieldnames(SolverTolerances))...)
    options=SolverOptions(α_Σ=0.1,max_scba=4,max_poisson=3,tolerances=tolerances,
        convergence=ConvergencePolicy(minimum_scba_iterations=4,
            required_consecutive_scba_passes=1,required_consecutive_poisson_passes=1,
            stagnation_window=0,stagnation_relative_improvement=0.0))
    algorithms=AlgorithmOptions(mixing=:anderson,anderson_history_depth=3,
        anderson_regularization=1e-10,anderson_damping=0.5)
    production=ProductionOptions(worker_count=1,algorithms=algorithms,checkpoint_every_scba=1,
        checkpoint_every_outer=1,progress_every_scba=0,progress_every_outer=0)
    execution=public_state_execution(production)
    states=ExactSolverState[]
    final=solve_production(problem;options,production_options=production,execution,
        checkpoint_sink=state->push!(states,state))
    bind=(;options,production_options=production,execution,
        domain_adaptation=DomainAdaptationPolicy(),kernel_options=nothing)
    # Rehashed payloads must still agree with their historical acceptance evidence.
    for state in filter(s->s.boundary in (:outer_accepted,:final),states)
        @test validate_exact_state(problem,state;bind...)===nothing
        rows=deepcopy(state.scba.history)
        old=last(rows)
        names=fieldnames(SCBAIteration)
        values=NamedTuple{names}(Tuple(getfield(old,name) for name in names))
        changed=merge(values,(;r_D=1e101))
        rows[end]=SCBAIteration((getfield(changed,name) for name in names)...)
        scba=state.scba
        wrong=SCBAResult(scba.green,scba.scattering,scba.embedding,scba.embedding_plus,
            scba.embedding_minus,rows,scba.converged,scba.status,scba.quality,
            scba.restart_contract,scba.mixer_state)
        @test_throws ArgumentError validate_exact_state(problem,
            replace_public_state(state;scba=wrong);bind...)
        p=state.progress
        wrongprogress=SolverProgress(p.completed_outer,p.current_outer,p.accepted_scba,
            p.outer_strict_streak,p.inner_strict_streak+1,p.diagnostic_streak)
        @test_throws ArgumentError validate_exact_state(problem,
            replace_public_state(state;progress=wrongprogress);bind...)
        if state.final_evidence!==nothing
            original=state.final_evidence
            receipt=FinalEvidenceSnapshot(report=original.report,
                original_converged=!original.original_converged,original_status=original.original_status,
                observations=original.observations)
            @test_throws ArgumentError validate_exact_state(problem,
                replace_public_state(state;final_evidence=receipt);bind...)
        end
    end
    # The same actual converged domain is terminal before its adaptive wrapper decision.
    terminal=only(filter(s->s.boundary===:final,states))
    policy=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=0,maximum_energy_nodes=31)
    kernels=ProductionKernelOptions()
    identity=exact_restart_identity(problem,options,production,execution;domain_adaptation=policy,kernel_options=kernels)
    adaptation=N._adaptation_snapshot(policy,kernels,identity,problem,algorithms,Dict{String,Any}[],
        0,:assess_measured_domain,nothing;elapsed_domain_seconds=0.0)
    domain=replace_public_state(terminal;identity,boundary=:domain_final,adaptation)
    adaptivebind=(;options,production_options=production,execution,domain_adaptation=policy,kernel_options=kernels)
    @test validate_exact_state(problem,domain;adaptivebind...)===nothing
    rows=deepcopy(domain.scba.history)
    names=fieldnames(SCBAIteration)
    values=NamedTuple{names}(Tuple(getfield(last(rows),name) for name in names))
    changed=merge(values,(;r_D=1e101))
    rows[end]=SCBAIteration((getfield(changed,name) for name in names)...)
    scba=domain.scba
    wrong=SCBAResult(scba.green,scba.scattering,scba.embedding,scba.embedding_plus,scba.embedding_minus,
        rows,scba.converged,scba.status,scba.quality,scba.restart_contract,scba.mixer_state)
    @test_throws ArgumentError validate_exact_state(problem,replace_public_state(domain;scba=wrong);adaptivebind...)
    useful=filter(state->state.boundary===:scba_accepted &&
        !state.continuation.decision.terminal && length(state.scba.mixer_state.states)>=2,states)
    @test !isempty(useful)
    if !isempty(useful)
        state=last(useful)
        c=state.continuation
        for populations in (c.next_previous_populations .+ 99.0,
                            [nextfloat(value) for value in c.next_previous_populations])
            wrong=SCBAContinuation(decision=c.decision,candidate=c.candidate,
                next_previous_current=c.next_previous_current,next_previous_populations=populations,
                seed_mu_scaled=c.seed_mu_scaled,seed_number_ratio=c.seed_number_ratio,
                witness_retention=c.witness_retention)
            @test_throws ArgumentError validate_exact_state(problem,
                replace_public_state(state;continuation=wrong);bind...)
        end
        old=state.scba
        for count in (0,length(old.mixer_state.states)-1)
            missing=SCBAMixerState(:anderson,old.mixer_state.states[1:count],old.mixer_state.residuals[1:count])
            wrong=SCBAResult(old.green,old.scattering,old.embedding,old.embedding_plus,old.embedding_minus,
                old.history,old.converged,old.status,old.quality,old.restart_contract,missing)
            @test_throws ArgumentError validate_exact_state(problem,
                replace_public_state(state;scba=wrong);bind...)
        end
    end
    selected=ExactSolverState[]
    isempty(useful) || push!(selected,last(useful))
    push!(selected,first(filter(s->s.boundary===:scba_accepted && s.continuation.decision.terminal,states)))
    push!(selected,first(filter(s->s.boundary===:outer_accepted,states)))
    for state in selected
        handle=restore_exact_state(problem,state;options,production_options=production,execution,
            domain_adaptation=DomainAdaptationPolicy(),kernel_options=nothing)
        rows=ScientificHistorySnapshot[]
        result=solve_production(problem;restart=handle,history_observer=row->push!(rows,row))
        @test result.status==final.status
        @test result.Uᴴ==final.Uᴴ
        @test result.scba.green.Gᴿ==final.scba.green.Gᴿ
        @test result.scba.green.Gˡ==final.scba.green.Gˡ
        for field in fieldnames(SCBAIteration)[1:17]
            @test isequal([getfield(row,field) for row in result.scba.history],
                [getfield(row,field) for row in final.scba.history])
        end
        for family in (:embedding,:embedding_plus,:embedding_minus), component in (:Σᴿ,:Σˡ,:Σᵍ)
            @test getfield(getfield(result.scba,family),component)==getfield(getfield(final.scba,family),component)
        end
        @test length(result.scba.mixer_state.states)==length(final.scba.mixer_state.states)
        @test all(row->!(row.action===:append && row.stage===:scba &&
            row.outer_iteration==state.progress.current_outer && row.inner_iteration<=state.progress.accepted_scba),rows)
    end
end

@testset "Actual guard reset and coefficient fallback preserve distinct mixer coverage" begin
    algorithms=AlgorithmOptions(mixing=:anderson,anderson_history_depth=3,
        anderson_regularization=0.0,anderson_damping=1.0)
    production=ProductionOptions(worker_count=1,algorithms=algorithms)
    workspace=N._AndersonWorkspace()
    state=[causal_family(4.0),causal_family(4.0)]
    N._anderson_mix!(workspace,state,[causal_family(1.0),causal_family(2.0)],algorithms,1.0,production)
    N._anderson_mix!(workspace,state,[causal_family(0.1),causal_family(1.5)],algorithms,1.0,production)
    # The accelerated first linewidth is -19/74; admissible convex fallback clears history.
    @test workspace.position.last_reset_reason===:causality_guard
    @test workspace.position.last_reset_mix==2
    @test workspace.position.applied_mixes==2
    @test isempty(workspace.position.retained_mix_indices) && isempty(workspace.states)
    mixer=SCBAMixerState(:anderson,deepcopy(workspace.states),deepcopy(workspace.residuals))
    @test N._validate_mixer_position(workspace.position,mixer,3,algorithms)===nothing
    resumed=N._AndersonWorkspace(deepcopy(mixer.states),deepcopy(mixer.residuals);position=workspace.position)
    other=deepcopy(state)
    candidate=[causal_family(0.2),causal_family(1.6)]
    N._anderson_mix!(workspace,state,candidate,algorithms,1.0,production)
    N._anderson_mix!(resumed,other,candidate,algorithms,1.0,production)
    @test workspace.position==resumed.position
    @test resumed.position.retained_mix_indices==(3,)
    @test all(j->all(name->getfield(state[j],name)==getfield(other[j],name),(:Σᴿ,:Σˡ,:Σᵍ)),1:2)
    fixed=N._AndersonWorkspace()
    for _ in 1:2
        N._anderson_mix!(fixed,[causal_family(2.0)],[causal_family(2.0)],algorithms,1.0,production)
    end
    @test fixed.position.last_coefficient_fallback===:coefficient_solve
    @test fixed.position.last_reset_reason===:none
    @test fixed.position.retained_mix_indices==(1,2)
    @test N._validate_mixer_position(fixed.position,
        SCBAMixerState(:anderson,fixed.states,fixed.residuals),3,algorithms)===nothing
    growing=N._AndersonWorkspace()
    N._anderson_mix!(growing,[causal_family(1.0)],[causal_family(2.0)],algorithms,0.25,production)
    N._anderson_mix!(growing,[causal_family(1.0)],[causal_family(5.0)],algorithms,0.25,production)
    @test growing.position.last_reset_reason===:residual_increase
    @test growing.position.retained_mix_indices==(2,)
    @test N._validate_mixer_position(growing.position,
        SCBAMixerState(:anderson,growing.states,growing.residuals),3,algorithms)===nothing
end
end
