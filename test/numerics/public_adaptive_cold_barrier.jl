module PublicAdaptiveColdBarrierTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")
const N=QCLNEGF.QCLNumerics
struct StopBeforeColdBuild <: Exception end
struct StopAfterAcceptedCheckpoint <: Exception end
function recorded_coverage(problem,coverage;measured=false)
    shared=(discretization_id=coverage.discretization_id,basis_id=coverage.basis_id,
        E_min_eV=ustrip(u"eV",problem.numerical.E_min),
        E_max_eV=ustrip(u"eV",problem.numerical.E_max),NE=problem.numerical.N_E,
        energy_min_scaled=coverage.energy_min,energy_max_scaled=coverage.energy_max,
        required_min_scaled=coverage.required_min,required_max_scaled=coverage.required_max,
        reason="explicit scalar replay fixture",seed_transfer="cold rebuild")
    evidence=measured ? DomainMeasuredEvidence(;shared...,hartree_id=coverage.hartree_id,
        tail_low=1.0,tail_high=0.0,edge_gamma=0.0,edge_spectral=0.0,scba_quality=:unresolved,
        seconds=0.0,status="explicit scalar replay fixture",strict_accepted=false,
        metrics=Dict{Symbol,Float64}(),current_A_per_m2=NaN,momentum_tail_unresolved=false) :
        DomainPreflightEvidence(;shared...)
    return DomainRevision(;revision=0,source=measured ? :measured_state : :preflight_estimate,evidence)
end

@testset "Recipe grids follow the owning ceil and IEEE unit roundtrip" begin
    f=public_state_fixture()
    policy=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=1,
        maximum_energy_nodes=31)
    kernels=ProductionKernelOptions()
    old=f.problem.numerical
    step=(old.E_max-old.E_min)/(old.N_E-1)
    numerical=NumericalParameters(; (name=>(name===:E_min ? old.E_min-step :
        name===:E_max ? old.E_max+step : name===:N_E ? 11 : getfield(old,name))
        for name in fieldnames(NumericalParameters))...)
    base=N._problem_recipe(f.problem,f.production.algorithms,kernels)
    recipe(n)=ProblemRebuildRecipe(; (name=>(name===:numerical ? n : getfield(base,name))
        for name in fieldnames(ProblemRebuildRecipe))...)
    # This literal 9 -> 11 recipe is a scalar geometry fixture, not an assertion
    # that a failed synthetic SCBA result could schedule a cold rebuild.
    @test_nowarn N._validate_problem_recipe(recipe(numerical),f.problem,f.production,kernels)
    for nodes in (3,9,15)
        changed=NumericalParameters(; (name=>(name===:N_E ? nodes : getfield(numerical,name))
            for name in fieldnames(NumericalParameters))...)
        @test_throws ArgumentError N._validate_problem_recipe(recipe(changed),f.problem,
            f.production,kernels)
    end
    shrink=NumericalParameters(; (name=>(name===:E_min ? old.E_min+step :
        name===:N_E ? 8 : getfield(old,name)) for name in fieldnames(NumericalParameters))...)
    @test_throws ArgumentError N._validate_problem_recipe(recipe(shrink),f.problem,f.production,kernels)
    changed_momentum=NumericalParameters(; (name=>(name===:N_k ? 3 : getfield(numerical,name))
        for name in fieldnames(NumericalParameters))...)
    @test_throws ArgumentError N._validate_problem_recipe(recipe(changed_momentum),f.problem,
        f.production,kernels)
    estimate=representation_coverage(f.problem)
    scaled_step=(estimate.energy_max-estimate.energy_min)/(old.N_E-1)
    values=Tuple(name===:required_min ? estimate.energy_min-scaled_step :
        name===:required_max ? estimate.energy_max+scaled_step :
        name===:covered ? false : getfield(estimate,name)
        for name in fieldnames(N.RepresentationCoverage))
    padded=N.RepresentationCoverage(values...)
    estimated=N._adapted_energy_parameters(f.problem,padded,policy)
    # Hand-checked IEEE case: scaled ceil gives 11; converting bounds to eV
    # before ceil gives 12 because the physical ratio is just above 10.
    @test estimated.N_E==11
    @test ceil(Int,ustrip(Unitful.NoUnits,(estimated.E_max-estimated.E_min)/step))+1==12
    @test_nowarn N._validate_problem_recipe(recipe(estimated),f.problem,f.production,kernels)
    grown=N._adapted_energy_parameters(f.problem,padded,policy;measured_failure=true)
    # Physical width 0.57 eV / 0.0475 eV rounds to 12.000000000000002.
    @test grown.N_E==14
    @test_nowarn N._validate_problem_recipe(recipe(grown),f.problem,f.production,kernels)
    @test N._recorded_domain_target(recorded_coverage(f.problem,padded),policy,f.problem,
        "preflight IEEE fixture")==((-1.0-0.47500000000000003)*0.1,(2.8000000000000003+0.47500000000000003)*0.1,11)
    @test N._recorded_domain_target(recorded_coverage(f.problem,padded;measured=true),policy,
        f.problem,"measured IEEE fixture")==(-0.195,0.375,14)
    for (measured,cap) in ((false,11),(true,14))
        at_cap=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=1,
            maximum_energy_nodes=cap)
        @test N._adapted_energy_parameters(f.problem,padded,at_cap;
            measured_failure=measured).N_E==cap
        below_cap=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=1,
            maximum_energy_nodes=cap-1)
        @test_throws ArgumentError N._adapted_energy_parameters(f.problem,padded,below_cap;
            measured_failure=measured)
    end
    pvalues=Tuple(name===:E_ref ? 41.0u"eV" : getfield(f.problem.physical,name)
        for name in fieldnames(PhysicalParameters))
    shifted_n=NumericalParameters(; (name=>(name===:E_min ? 40.9u"eV" :
        name===:E_max ? 41.28u"eV" : getfield(old,name))
        for name in fieldnames(NumericalParameters))...)
    shifted=build_problem(physical=PhysicalParameters(pvalues...),numerical=shifted_n,
        scales=ScaleSystem(E₀=0.07u"eV"),scattering=f.problem.scattering,
        physical_models=f.problem.models)
    coverage=representation_coverage(shifted)
    shifted_step=(coverage.energy_max-coverage.energy_min)/(old.N_E-1)
    values=Tuple(name===:required_min ? coverage.energy_min-0.125shifted_step :
        name===:required_max ? coverage.energy_max+0.25shifted_step :
        name===:covered ? false : getfield(coverage,name)
        for name in fieldnames(N.RepresentationCoverage))
    shifted_next=N._adapted_energy_parameters(shifted,N.RepresentationCoverage(values...),policy)
    @test shifted_next.N_E==10 # 8 original intervals plus 0.375 requires 9 intervals.
    row=recorded_coverage(shifted,N.RepresentationCoverage(values...))
    @test isequal(N._recorded_domain_target(row,policy,shifted,"nonzero reference fixture"),
        (ustrip(u"eV",shifted_next.E_min),ustrip(u"eV",shifted_next.E_max),10))
    shifted_base=N._problem_recipe(shifted,f.production.algorithms,kernels)
    target=ProblemRebuildRecipe(; (name=>(name===:numerical ? shifted_next :
        getfield(shifted_base,name)) for name in fieldnames(ProblemRebuildRecipe))...)
    @test_nowarn N._validate_problem_recipe(target,shifted,f.production,kernels)
    coarser=NumericalParameters(; (name=>(name===:N_E ? 9 : getfield(shifted_next,name))
        for name in fieldnames(NumericalParameters))...)
    bad=ProblemRebuildRecipe(; (name=>(name===:numerical ? coarser : getfield(target,name))
        for name in fieldnames(ProblemRebuildRecipe))...)
    @test_throws ArgumentError N._validate_problem_recipe(bad,shifted,f.production,kernels)
end

@testset "Historical preflight grids and typed policy agree with the current domain" begin
    f=public_state_fixture()
    policy=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=1,
        maximum_energy_nodes=31)
    kernels=ProductionKernelOptions()
    identity=exact_restart_identity(f.problem,f.options,f.production,f.execution;
        domain_adaptation=policy,kernel_options=kernels)
    current=N._problem_recipe(f.problem,f.production.algorithms,kernels)
    estimate=representation_coverage(f.problem)
    preflight=DomainPreflightEvidence(discretization_id="synthetic previous grid",
        basis_id=estimate.basis_id,E_min_eV=0.0,E_max_eV=0.19000000000000003,NE=5,
        energy_min_scaled=0.0,energy_max_scaled=1.9000000000000001,
        required_min_scaled=-1.0,required_max_scaled=2.8000000000000003,
        reason="synthetic preflight proposal",seed_transfer="none; no state solved")
    function after_preflight(row;changes...)
        adaptation=DomainAdaptationState(;merge((;policy,kernel_options=kernels,expansions=1,
            elapsed_domain_seconds=0.0,current_problem=identity.problem,current_problem_recipe=current,
            revisions=(DomainRevision(revision=0,source=:preflight_estimate,evidence=row),),
            next_action=:assess_measured_domain,decision_id=1,decision_committed=false,
            next_problem_recipe=nothing),(;changes...))...)
        return replace_public_state(f.state;identity,adaptation)
    end
    bind=(options=f.options,production_options=f.production,execution=f.execution,
        domain_adaptation=policy,kernel_options=kernels)
    @test validate_exact_state(f.problem,after_preflight(preflight);bind...)===nothing
    for (field,value) in ((:NE,3),(:NE,2),(:E_min_eV,NaN),(:E_max_eV,-0.1),
        (:required_max_scaled,3.8),(:required_min_scaled,Inf))
        row=DomainPreflightEvidence(; (name=>(name===field ? value : getfield(preflight,name))
            for name in fieldnames(DomainPreflightEvidence))...)
        @test_throws ArgumentError validate_exact_state(f.problem,after_preflight(row);bind...)
    end
    covered=DomainPreflightEvidence(; (name=>(name===:required_min_scaled ? 0.1 :
        name===:required_max_scaled ? 0.15 : getfield(preflight,name))
        for name in fieldnames(DomainPreflightEvidence))...)
    @test_throws ArgumentError validate_exact_state(f.problem,after_preflight(covered);bind...)
    for cap in (9,8)
        capped=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=1,
            maximum_energy_nodes=cap)
        capped_identity=exact_restart_identity(f.problem,f.options,f.production,f.execution;
            domain_adaptation=capped,kernel_options=kernels)
        checkpoint=replace_public_state(after_preflight(preflight;policy=capped);
            identity=capped_identity)
        capped_bind=merge(bind,(;domain_adaptation=capped))
        if cap==9
            @test validate_exact_state(f.problem,checkpoint;capped_bind...)===nothing
        else
            @test_throws ArgumentError validate_exact_state(f.problem,checkpoint;capped_bind...)
        end
    end
    wrong_policy=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=1,
        maximum_energy_nodes=31,growth_fraction=nextfloat(policy.growth_fraction))
    @test_throws ArgumentError validate_exact_state(f.problem,
        after_preflight(preflight;policy=wrong_policy);bind...)
    @test_throws ArgumentError validate_exact_state(f.problem,
        after_preflight(preflight;revisions=());bind...)
    @test_throws ArgumentError validate_exact_state(f.problem,
        after_preflight(preflight;decision_committed=true);bind...)
    mismatch=DomainRevision(0,:measured_state,preflight)
    @test_throws ArgumentError validate_exact_state(f.problem,
        after_preflight(preflight;revisions=(mismatch,));bind...)
end

@testset "Actual owner cold decision validates and resumes exactly one new domain" begin
    f=public_state_fixture()
    # Explicit engineering inputs keep both domains small. The tiny tail threshold
    # selects the adaptive branch; no physical convergence/accuracy is claimed.
    numerical=NumericalParameters(; (name=>(name===:k_max ? 0.1u"nm^-1" :
        name===:E_min ? -1.0u"eV" : name===:E_max ? 1.0u"eV" :
        getfield(f.problem.numerical,name)) for name in fieldnames(NumericalParameters))...)
    problem=build_problem(physical=reference_parameters(f_ion=1e-8),numerical=numerical,
        scales=f.problem.scales,scattering=f.problem.scattering,physical_models=f.problem.models)
    covered=representation_coverage(problem).covered
    @test covered
    covered || error("engineering fixture must cover its Hamiltonian before the bounded solve")
    options=SolverOptions(max_scba=1,max_poisson=1)
    production=f.production
    execution=f.execution
    policy=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=1,
        maximum_energy_nodes=13,tail_threshold=1e-100)
    kernels=ProductionKernelOptions()
    saved=Ref{Union{Nothing,ExactSolverState}}(nothing)
    original_history=ScientificHistorySnapshot[]
    @test_throws StopBeforeColdBuild solve_adaptive_production(problem;options,
        production_options=production,execution,domain_adaptation=policy,kernel_options=kernels,
        history_observer=row->push!(original_history,row),
        checkpoint_sink=state->begin
            if state.boundary===:domain_final && state.adaptation.next_action===:cold_rebuild
                saved[]=state
                throw(StopBeforeColdBuild())
            end
        end)
    @test saved[]!==nothing
    if saved[]!==nothing
        state=saved[]
        adaptation=state.adaptation
        target=adaptation.next_problem_recipe
        revision=only(adaptation.revisions)
        bind=(;options,production_options=production,execution,domain_adaptation=policy,
            kernel_options=kernels)
        @test state.scba.quality!==:invalid
        @test target.numerical.N_E==13 # 3 eV widened window / 0.25 eV step = 12 intervals.
        @test validate_exact_state(problem,state;bind...)===nothing
        @test !isempty(original_history) && all(row->row.domain_revision==0 &&
            row.discretization_id==revision.evidence.discretization_id &&
            row.basis_id==revision.evidence.basis_id,original_history)
        function with_adaptation(checkpoint;changes...)
            changed=DomainAdaptationState(; (name=>(haskey(changes,name) ? changes[name] :
                getfield(checkpoint.adaptation,name)) for name in fieldnames(DomainAdaptationState))...)
            return replace_public_state(checkpoint;adaptation=changed)
        end
        for nodes in (3,9,19)
            changed=NumericalParameters(; (name=>(name===:N_E ? nodes : getfield(target.numerical,name))
                for name in fieldnames(NumericalParameters))...)
            recipe=ProblemRebuildRecipe(; (name=>(name===:numerical ? changed : getfield(target,name))
                for name in fieldnames(ProblemRebuildRecipe))...)
            @test_throws ArgumentError validate_exact_state(problem,
                with_adaptation(state;next_problem_recipe=recipe);bind...)
        end
        # Same preserved 0.25 eV step but only half of the declared growth.
        # Geometry alone cannot distinguish this from the owning proposal.
        undergrown=NumericalParameters(; (name=>(name===:E_min ? -1.25u"eV" :
            name===:E_max ? 1.25u"eV" : name===:N_E ? 11 : getfield(target.numerical,name))
            for name in fieldnames(NumericalParameters))...)
        undergrown_recipe=ProblemRebuildRecipe(; (name=>(name===:numerical ? undergrown :
            getfield(target,name)) for name in fieldnames(ProblemRebuildRecipe))...)
        @test_throws ArgumentError validate_exact_state(problem,
            with_adaptation(state;next_problem_recipe=undergrown_recipe);bind...)
        premature=replace_public_state(with_adaptation(state;next_action=:finished,
            next_problem_recipe=nothing);boundary=:final)
        @test_throws ArgumentError validate_exact_state(problem,premature;bind...)
        for (field,value) in ((:NE,3),(:E_min_eV,nextfloat(revision.evidence.E_min_eV)),
            (:discretization_id,"wrong current grid"),(:basis_id,"wrong current basis"),
            (:hartree_id,"wrong published Hartree"),(:strict_accepted,!revision.evidence.strict_accepted),
            (:energy_min_scaled,nextfloat(revision.evidence.energy_min_scaled)),
            (:required_max_scaled,Inf),(:scba_quality,:unknown))
            row=DomainMeasuredEvidence(; (name=>(name===field ? value : getfield(revision.evidence,name))
                for name in fieldnames(DomainMeasuredEvidence))...)
            @test_throws ArgumentError validate_exact_state(problem,
                with_adaptation(state;revisions=(DomainRevision(revision=0,source=:measured_state,evidence=row),));bind...)
        end
        # Same coherent terminal history/receipt failure, but invalid SCBA cannot
        # schedule a rebuild: the owning loop always finishes that domain.
        failed=public_state_failed_scba(state.scba)
        status=:scba_invalid_candidate
        evidence=FinalEvidenceSnapshot(report=state.final_evidence.report,original_converged=false,
            original_status=status,observations=state.final_evidence.observations)
        row=DomainMeasuredEvidence(; (name=>(name===:status ? String(status) :
            name===:scba_quality ? :invalid :
            getfield(revision.evidence,name)) for name in fieldnames(DomainMeasuredEvidence))...)
        invalid=replace_public_state(state;scba=failed,status,final_evidence=evidence)
        invalid=with_adaptation(invalid;
            revisions=(DomainRevision(revision=0,source=:measured_state,evidence=row),))
        @test_throws ArgumentError validate_exact_state(problem,invalid;bind...)
        # The same failed terminal state is legitimately finished by the actual
        # invalid-quality branch even with budget left. Failed IEEE tails remain raw.
        failed_row=DomainMeasuredEvidence(; (name=>(name===:tail_low ? NaN :
            name===:tail_high ? Inf : getfield(row,name)) for name in fieldnames(DomainMeasuredEvidence))...)
        finished_invalid=replace_public_state(with_adaptation(invalid;next_action=:finished,
            next_problem_recipe=nothing,
            revisions=(DomainRevision(revision=0,source=:measured_state,evidence=failed_row),));boundary=:final)
        @test validate_exact_state(problem,finished_invalid;bind...)===nothing
        @test isnan(only(finished_invalid.adaptation.revisions).evidence.tail_low)
        @test only(finished_invalid.adaptation.revisions).evidence.tail_high===Inf
        budget=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=0,
            maximum_energy_nodes=13,tail_threshold=policy.tail_threshold)
        budget_identity=exact_restart_identity(problem,options,production,execution;
            domain_adaptation=budget,kernel_options=kernels)
        finished_budget=replace_public_state(with_adaptation(premature;policy=budget);
            identity=budget_identity)
        @test validate_exact_state(problem,finished_budget;
            merge(bind,(;domain_adaptation=budget))...)===nothing

        handle=restore_exact_state(problem,state;bind...)
        next=Ref{Union{Nothing,ExactSolverState}}(nothing)
        history=ScientificHistorySnapshot[]
        @test_throws StopAfterAcceptedCheckpoint solve_adaptive_production(problem;restart=handle,
            checkpoint_sink=s->begin next[]=s; throw(StopAfterAcceptedCheckpoint()); end,
            history_observer=row->push!(history,row))
        @test next[]!==nothing
        if next[]!==nothing
            @test next[].boundary===:scba_accepted
            @test next[].progress.accepted_scba==1
            @test next[].progress.current_outer==1
            @test next[].adaptation.expansions==1
            @test length(next[].adaptation.revisions)==1
            @test next[].adaptation.revisions[1].revision==0
            @test next[].adaptation.current_problem_recipe.numerical.N_E==13
            @test next[].adaptation.current_problem_recipe.numerical.E_min===target.numerical.E_min
            @test next[].adaptation.kernel_options===kernels
            @test all(row->row.stage===:scba && row.inner_iteration==1,history)
            rebuilt=N._rebuild_from_recipe(target) # Same small recipe, no additional solve.
            coverage=representation_coverage(rebuilt)
            @test !isempty(history) && all(row->row.domain_revision==1 &&
                row.discretization_id==coverage.discretization_id && row.basis_id==coverage.basis_id,history)
            @test validate_exact_state(rebuilt,next[];bind...)===nothing
            for (field,value) in ((:NE,3),(:E_max_eV,2.0),(:scba_quality,:invalid))
                row=DomainMeasuredEvidence(; (name=>(name===field ? value :
                    getfield(revision.evidence,name)) for name in fieldnames(DomainMeasuredEvidence))...)
                changed=with_adaptation(next[];
                    revisions=(DomainRevision(revision=0,source=:measured_state,evidence=row),))
                @test_throws ArgumentError validate_exact_state(rebuilt,changed;bind...)
            end
        end
        # No checkpoint factory exists in this invocation. The actual resumed
        # domain must still be bound to every history row, with revision 1.
        history_only=ScientificHistorySnapshot[]
        result=solve_adaptive_production(problem;restart=handle,
            history_observer=row->push!(history_only,row))
        coverage=representation_coverage(result.problem)
        @test !isempty(history_only) && all(row->row.domain_revision==1 &&
            row.discretization_id==coverage.discretization_id && row.basis_id==coverage.basis_id,history_only)
        @test all(row->row.discretization_id===first(history_only).discretization_id &&
            row.basis_id===first(history_only).basis_id,history_only)
        @test last(result.observables[:domain_adaptation]["revisions"])["revision"]==1
    end
end
end
