isdefined(@__MODULE__,:FinalObservationContractTests) || include("../domain/final_observation_contract.jl")
module FinalObservationSnapshotTests
include("../support/common.jl")
using ..FinalObservationContractTests: final_fixture,replace_solution,take_snapshot,caught,replace_public_state,public_state_execution
const N=QCLNEGF.QCLNumerics

function reseal(f;solution=f.solution)
    id=exact_restart_identity(solution.problem,f.options,f.production,f.execution;
        domain_adaptation=f.policy,kernel_options=nothing)
    state=replace_public_state(f.state;identity=id,published_U=solution.Uᴴ,published_n=solution.n,
        scba=solution.scba,outer_history=solution.outer_history,final_evidence=N._final_evidence(solution))
    return (;f...,problem=solution.problem,solution,state)
end
function replace_problem(p;changes...)
    names=fieldnames(typeof(p))
    args=merge(NamedTuple{names}(Tuple(getfield(p,n) for n in names)),(;changes...))
    NEGFProblem((getfield(args,n) for n in names)...)
end

@testset "Stored full coherences, units and display gauge have literal oracles" begin
    f=final_fixture(basis_states=2);p=f.problem;g=p.grids;b=p.basis;s=p.scales
    # Manufacture a failed-state observation packet, not a physical accuracy fixture.
    g.wᴱ .= [0.13,0.17,0.19,0.23,0.29,0.31,0.37,0.41,0.43]
    g.wᵏ .= [0.07,0.11]
    g.wˣ .*= [1+i/100 for i in eachindex(g.wˣ)]
    for k in eachindex(g.κ),e in eachindex(g.ε)
        matrix=ComplexF64[0.5+e/10 (0.2+0.1im)*k;(0.2-0.1im)*k -0.1]
        f.solution.scba.green.Gˡ[e,k,:,:].=im.*matrix
        f.solution.scba.green.A[e,k,:,:].=ComplexF64[2+e/10 0.3im;-0.3im 1+k/10]
    end
    f=reseal(f)
    before=N._exact_digest((p,f.solution,f.state))
    snap=take_snapshot(f);v=snap.views;c=snap.context
    @test snap.summary.scientific_accepted===false
    @test snap.summary.quality===:invalid
    @test isinf(snap.metrics.summary_metrics[:keldysh]) || isnan(snap.metrics.summary_metrics[:keldysh])
    @test c.dimensions==(25,9,2,2)
    @test c.domain_revision==0
    @test c.inputs.g_s isa Int
    @test c.state_index==Int64[0,1]
    @test c.wE==g.wᴱ && c.wk==g.wᵏ && c.wX==g.wˣ
    @test c.boundary_condition===:stationary_field_periodic_embedding
    @test c.energy_eV==g.ε .* s.E₀_eV .+ ustrip(u"eV",p.physical.E_ref)
    @test c.z_nm==g.x .* s.L₀_m .* 1e9
    @test c.k_per_nm==g.κ ./ s.L₀_m .* 1e-9
    # Independent scalar sums, no observable getter in expected path.
    expectedP=zeros(ComplexF64,2,2)
    expectedEz=zeros(9,25)
    for a in 1:2,bb in 1:2,e in 1:9,k in 1:2
        expectedP[a,bb]+=(-im*p.physical.g_s/(2π*s.L₀_m^2))*g.wᴱ[e]*g.wᵏ[k]*
            f.solution.scba.green.Gˡ[e,k,a,bb]
    end
    for z in 1:25,e in 1:9,k in 1:2,a in 1:2,bb in 1:2
        expectedEz[e,z]+=p.physical.g_s/(2π*s.E₀_eV*s.L₀_m^3)*g.wᵏ[k]*
            real(-im*b.χ[z,a]*f.solution.scba.green.Gˡ[e,k,a,bb]*conj(b.χ[z,bb]))
    end
    @test v.sheet_density_matrix_per_m2≈expectedP rtol=5e-15
    @test v.spatial_energy_density_per_eV_m3≈expectedEz rtol=5e-15
    @test imag(v.sheet_density_matrix_per_m2[1,2])!=0
    for e in 1:9,k in 1:2
        @test v.spectral_trace_per_eV[e,k]≈(3+e/10+k/10)/s.E₀_eV rtol=5e-15
        @test v.occupied_trace_per_eV[e,k]≈(0.4+e/10)/s.E₀_eV rtol=5e-15
    end
    # Literal legacy display procedure, separate from the getter and from public spectrum gate.
    h0=Matrix(view(project_hamiltonians(p,f.solution.Uᴴ),1,:,:))
    eigenvalues=eigvals(Hermitian(h0))
    expected_levels=[real(tr(h0))/2-sqrt(((real(h0[1,1])-real(h0[2,2]))/2)^2+abs2(h0[1,2])),
        real(tr(h0))/2+sqrt(((real(h0[1,1])-real(h0[2,2]))/2)^2+abs2(h0[1,2]))]
    @test v.effective_levels_eV≈expected_levels .* s.E₀_eV .+ c.inputs.E_ref_eV
    @test eigenvalues≈expected_levels
    @test h0*v.effective_transform≈v.effective_transform*Diagonal(expected_levels)
    @test v.effective_density_matrix_per_m2≈v.effective_transform'*expectedP*v.effective_transform
    @test v.effective_wavefunctions_per_sqrt_nm≈b.χ*v.effective_transform/sqrt(s.L₀_m*1e9)
    for column in eachcol(v.effective_wavefunctions_per_sqrt_nm)
        pivot=argmax(abs.(column))
        @test abs(imag(column[pivot]))<1e-12
        @test real(column[pivot])>=0
    end
    @test v.potential_total_eV==v.potential_structure_eV.+v.potential_external_eV.+v.potential_hartree_eV
    @test v.density_per_m3==f.solution.n./s.L₀_m^3
    @test v.sheet_density_per_m2≈real(tr(expectedP))
    # Literal four-array audit integrals, including a negative coherence mode.
    integral=zeros(ComplexF64,2,2,2);cumulative=zeros(9,2)
    for k in 1:2,e in 1:9
        for a in 1:2,bb in 1:2
            integral[k,a,bb]+=g.wᴱ[e]/(2π)*f.solution.scba.green.A[e,k,a,bb]
        end
        cumulative[e,k]=real(integral[k,1,1]+integral[k,2,2])
    end
    @test snap.positivity.scan_order===:momentum_energy
    @test snap.positivity.integrals.spectral_integral≈integral
    @test snap.positivity.integrals.spectral_trace_cumulative≈cumulative
    @test snap.positivity.integrals.total_blocks_per_kind==18
    @test N._exact_digest((p,f.solution,f.state))==before
    # Outputs detach every owned array, and the original public audit metadata.
    snap.context.wE[1]=99;snap.views.sheet_density_matrix_per_m2[1]=123
    snap.positivity.integrals.spectral_integral[1]=456
    @test N._exact_digest((p,f.solution,f.state))==before
end

@testset "Failed IEEE state remains raw; unmeasurable views have exact issues" begin
    f=final_fixture()
    f.solution.Uᴴ[1]=-0.0
    payload=reinterpret(Float64,0x7ff8000000000042)
    f.solution.scba.green.Gˡ[1]=complex(payload,0.0)
    f=reseal(f)
    snap=take_snapshot(f)
    @test snap.views.spatial_energy_density_per_eV_m3===nothing
    @test only(filter(i->i.target===:spatial_energy_density_per_eV_m3,snap.issues)).reason===:nonfinite_green
    @test isnan(snap.views.occupied_trace_per_eV[1])
    @test any(w->w.status===:error,snap.positivity.witnesses)
    wrong=deepcopy(f.solution);wrong.Uᴴ[1]=0.0
    @test caught(()->take_snapshot((;f...,solution=wrong))) isa ArgumentError
    wrong=deepcopy(f.solution);wrong.scba.green.Gˡ[1]=complex(reinterpret(Float64,0x7ff8000000000043),0.0)
    @test caught(()->take_snapshot((;f...,solution=wrong))) isa ArgumentError
    @test reinterpret(UInt64,real(f.state.scba.green.Gˡ[1]))==0x7ff8000000000042
    # The nonparabolic projection owner explicitly refuses nonfinite U; snapshot must
    # preserve a valid failed terminal rather than require an eigensystem to save it.
    p=replace_problem(f.problem;models=PhysicalModelOptions(dispersion=:kane_inplane,
        nonparabolicity_per_eV=0.02))
    sol=replace_solution(f.solution;problem=p);sol.Uᴴ[2]=NaN
    bad_projection=reseal(f;solution=sol)
    captured=take_snapshot(bad_projection)
    @test captured.views.effective_transform===nothing
    @test any(i->i.target===:effective_levels_eV,captured.issues)
end

@testset "Completed-parent receipts remain separate without evaluating a map" begin
    f=final_fixture()
    report=ConvergenceReport(false,Dict{Symbol,Float64}(),["manufactured unmeasurable revalidation"])
    ctx=AcceptanceContext(backend=:not_executed,algorithms=f.production.algorithms,
        hilbert_roundoff_measured=nothing,map_kernels_shared_with_solver=nothing,
        final_state_alignment=:restored_aligned)
    reval=RevalidationReceipt(status=:error,reason=:nonfinite_terminal_state,detail="not executed",
        context=ctx,report=report)
    oracle=IndependentOracleReceipt(status=:unsupported,reason=:independent_map_not_qualified,
        detail="constructor packet",context=nothing,report=nothing)
    obs=deepcopy(f.solution.observables)
    obs[:restoration_evidence]=CompletedStateEvidence(original=f.state.final_evidence,
        revalidation=reval,independent_oracle=oracle)
    obs[:acceptance_measurement_context]=Dict("backend"=>"not_executed",
        "hilbert_roundoff_measured"=>false,"map_kernels_shared_with_solver"=>nothing,
        "final_state_alignment"=>"restored_aligned")
    obs[:quality]=:invalid;obs[:termination_reason]="validation_failed"
    sol=replace_solution(f.solution;observables=obs,report,status=:validation_failed)
    snap=take_snapshot((;f...,solution=sol);source_mode=:completed_restore)
    @test snap.revalidation.status===:error
    @test snap.independent_oracle.status===:unsupported
    @test isnan(f.state.final_evidence.report.metrics[:r_K])
    @test isempty(snap.revalidation.report.metrics)
    @test !snap.summary.stationary_candidate_accepted
    @test !snap.summary.scientific_accepted
    @test caught(()->take_snapshot((;f...,solution=sol))) isa ArgumentError
    @test caught(()->take_snapshot(f;source_mode=:completed_restore)) isa ArgumentError
    wrong=deepcopy(sol);wrong.report.messages[1]="other solution"
    @test caught(()->take_snapshot((;f...,solution=wrong);source_mode=:completed_restore)) isa ArgumentError
    fabricated=ConvergenceReport(true,Dict{Symbol,Float64}(),["forged strict success"])
    false_observations=deepcopy(obs)
    false_observations[:restoration_evidence]=CompletedStateEvidence(original=f.state.final_evidence,
        revalidation=RevalidationReceipt(status=:error,reason=:nonfinite_terminal_state,
            detail="not executed",context=ctx,report=fabricated),independent_oracle=oracle)
    false_success=replace_solution(sol;report=fabricated,observables=false_observations,
        converged=true,status=:converged)
    @test caught(()->take_snapshot((;f...,solution=false_success);source_mode=:completed_restore)) isa ArgumentError
    snap.revalidation.report.messages[1]="recipient owns this"
    @test reval.report.messages[1]=="manufactured unmeasurable revalidation"
end

function failed_solution_on_problem(f,p)
    n=p.numerical;shape=(n.N_E,n.N_k,n.N_b,n.N_b)
    green=GreenState(fill(1.0-1im,shape),fill(0.5im,shape),fill(-1.5im,shape),
        fill(2.0+0im,shape),ones(n.N_E,n.N_k),ones(n.N_E,n.N_k))
    family()=SelfEnergyFamily(zeros(ComplexF64,shape),zeros(ComplexF64,shape),zeros(ComplexF64,shape))
    old=f.solution.scba
    scba=SCBAResult(green,Dict(name=>family() for name in p.kernels.enabled),family(),family(),family(),
        old.history,old.converged,old.status,old.quality,old.restart_contract,old.mixer_state)
    replace_solution(f.solution;problem=p,scba,Uᴴ=zeros(n.N_z),n=ones(n.N_z))
end

@testset "Actual zero IFR and EE families preserve construction flags and distinct model context" begin
    f=final_fixture()
    physical=PhysicalParameters((name===:Δᴵᶠᴿ ? 0.0u"nm" : getfield(f.problem.physical,name)
        for name in fieldnames(PhysicalParameters))...)
    ee=ElectronElectronOptions(mode=:sppa_single_q,screening_dimension=2,
        transfer_wavenumber=0.1u"nm^-1",electron_temperature=100u"K",
        effective_mass_ratio=0.067,carrier_density=1e10u"cm^-2")
    flags=ScatteringOptions(LO=false,acoustic=false,impurity=false,IFR=true,alloy=false)
    p=build_problem_production(;physical,numerical=f.problem.numerical,scattering=flags,
        physical_models=PhysicalModelOptions(electron_electron=ee),production_options=f.production).problem
    @test p.scattering.IFR
    @test !(:IFR in p.kernels.enabled)
    @test p.kernels.enabled==[:electron_electron]
    solution=failed_solution_on_problem(f,p)
    id=exact_restart_identity(p,f.options,f.production,f.execution;
        domain_adaptation=f.policy,kernel_options=nothing)
    state=replace_public_state(f.state;identity=id,published_U=solution.Uᴴ,published_n=solution.n,
        scba=solution.scba,family_order=(:electron_electron,:embedding_plus,:embedding_minus),
        final_evidence=N._final_evidence(solution))
    packet=(;f...,problem=p,solution,state)
    snap=take_snapshot(packet)
    @test snap.context.problem_scattering.IFR
    @test snap.context.enabled_scattering==()
    @test snap.context.selfenergy_families==(:electron_electron,)
    @test snap.context.physical_models.electron_electron.mode===:sppa_single_q
    @test all(name->getfield(snap.context.problem_scattering,name)==getfield(flags,name),fieldnames(ScatteringOptions))
    # Same-problem substituted per-mechanism tensor must be rejected, not hidden by an empty five-family list.
    wrong=deepcopy(solution);wrong.scba.scattering[:electron_electron].Σˡ[1]=im
    @test caught(()->take_snapshot((;packet...,solution=wrong))) isa ArgumentError
end

@testset "Adaptive final context binds the actual rebuilt problem; domain final is insufficient" begin
    f=final_fixture();p=f.problem
    policy=DomainAdaptationPolicy(mode=:expand_energy_window,maximum_expansions=1,maximum_energy_nodes=4097)
    kernels=ProductionKernelOptions()
    bind=(;f.binding...,domain_adaptation=policy,kernel_options=kernels)
    identity=exact_restart_identity(p,f.options,f.production,f.execution;
        domain_adaptation=policy,kernel_options=kernels)
    recipe(p)=ProblemRebuildRecipe(physical=p.physical,numerical=p.numerical,scales=p.scales,
        scattering=p.scattering,models=p.models,algorithms=f.production.algorithms,kernel_options=kernels)
    pending=DomainAdaptationState(policy=policy,kernel_options=kernels,expansions=0,
        elapsed_domain_seconds=0.0,current_problem=identity.problem,current_problem_recipe=recipe(p),
        revisions=(),next_action=:assess_measured_domain,decision_id=0,decision_committed=false,next_problem_recipe=nothing)
    domain=replace_public_state(f.state;identity,boundary=:domain_final,adaptation=pending)
    @test validate_exact_state(p,domain;bind...)===nothing
    @test caught(()->take_snapshot((;f...,state=domain,binding=bind))) isa ArgumentError
    # Explicit bounded coverage/postprocessing constructs receipts. This is not a new SCBA/Poisson solve.
    old=representation_coverage(p;hartree=f.solution.Uᴴ)
    @test !old.covered
    numerical=expand_energy_window(p,old;maximum_nodes=4097)
    current=build_problem(;physical=p.physical,numerical,scattering=p.scattering,scales=p.scales)
    solution=failed_solution_on_problem(f,current)
    coverage=representation_coverage(current;hartree=solution.Uᴴ)
    preflight=DomainPreflightEvidence(discretization_id=old.discretization_id,basis_id=old.basis_id,
        E_min_eV=ustrip(u"eV",p.numerical.E_min),E_max_eV=ustrip(u"eV",p.numerical.E_max),NE=p.numerical.N_E,
        energy_min_scaled=old.energy_min,energy_max_scaled=old.energy_max,
        required_min_scaled=old.required_min,required_max_scaled=old.required_max,
        reason="constructor-only source preflight",seed_transfer="none; no state solved")
    measured=DomainMeasuredEvidence(discretization_id=coverage.discretization_id,basis_id=coverage.basis_id,
        hartree_id=coverage.hartree_id,E_min_eV=ustrip(u"eV",current.numerical.E_min),
        E_max_eV=ustrip(u"eV",current.numerical.E_max),NE=current.numerical.N_E,
        energy_min_scaled=coverage.energy_min,energy_max_scaled=coverage.energy_max,
        required_min_scaled=coverage.required_min,required_max_scaled=coverage.required_max,
        tail_low=NaN,tail_high=Inf,edge_gamma=NaN,edge_spectral=Inf,scba_quality=:invalid,
        seconds=0.0,status=String(solution.status),strict_accepted=false,metrics=solution.report.metrics,
        current_A_per_m2=NaN,momentum_tail_unresolved=true,seed_transfer="cold rebuild",reason="failed raw evidence")
    identity=exact_restart_identity(current,f.options,f.production,f.execution;
        domain_adaptation=policy,kernel_options=kernels)
    adaptation=DomainAdaptationState(policy=policy,kernel_options=kernels,expansions=1,
        elapsed_domain_seconds=0.0,current_problem=identity.problem,current_problem_recipe=recipe(current),
        revisions=(DomainRevision(revision=0,source=:preflight_estimate,evidence=preflight),
            DomainRevision(revision=1,source=:measured_state,evidence=measured)),
        next_action=:finished,decision_id=2,decision_committed=true,next_problem_recipe=nothing)
    state=replace_public_state(f.state;identity,published_U=solution.Uᴴ,published_n=solution.n,
        scba=solution.scba,adaptation,final_evidence=N._final_evidence(solution))
    @test validate_exact_state(current,state;bind...)===nothing
    snap=take_snapshot((;f...,problem=current,solution,state,binding=bind))
    @test snap.context.domain_revision==1
    @test snap.context.dimensions[2]==current.numerical.N_E>p.numerical.N_E
    @test snap.context.discretization_id==coverage.discretization_id!=old.discretization_id
    @test snap.context.energy_scaled==current.grids.ε
    @test caught(()->final_observation_snapshot(replace_solution(solution;problem=p),state;bind...,
        source_mode=:produced,maximum_bytes=10^8,maximum_metadata_bytes=10^7)) isa Exception
end

@testset "Nested audit and all input-owned outputs are detached" begin
    f=final_fixture()
    f.solution.observables[:final_fixed_point_audit]=Dict("diagnostic"=>Float64[1,2])
    f=reseal(f)
    snap=take_snapshot(f)
    old=snap.context.wE[1]
    f.problem.grids.wᴱ[1]=old+1
    @test snap.context.wE[1]==old
    @test snap.assessment["quadrature"]["energy"]["weights"][1]==old
    snap.assessment["final_fixed_point_audit"]["diagnostic"][1]=99
    @test f.solution.observables[:final_fixed_point_audit]["diagnostic"][1]==1
    snap.metrics.summary_metrics[:keldysh]=123
    @test isnan(f.solution.report.metrics[:r_K])
end

@testset "Integrated PSD preserves its first-on-tie order and all four arrays" begin
    f=final_fixture()
    f.solution.scba.green.A .=1
    f.solution.scba.green.A[2,1,1,1]=-2
    f.solution.scba.green.A[1,2,1,1]=-2
    f=reseal(f)
    snap=take_snapshot(f);d=snap.positivity
    @test (d.witnesses[1].energy_index,d.witnesses[1].momentum_index)==(2,1)
    for k in 1:2
        expected=sum(f.problem.grids.wᴱ[e]/(2π)*real(f.solution.scba.green.A[e,k,1,1]) for e in 1:9)
        @test d.integrals.spectral_integral[k,1,1]≈expected
        @test d.integrals.spectral_integral_eigenvalues[k,1]≈expected
        negative=2*f.problem.grids.wᴱ[k==1 ? 2 : 1]/(2π)
        @test d.integrals.negative_energy_integral_scaled[1,k]≈negative
    end
end

# Report/history classification control only. These zeros are manufactured
# evidence, not a map measurement or a claim that the supplied G solves physics.
function classified_fixture(;approximate=false)
    f=final_fixture()
    policy=ConvergencePolicy(mode=:adaptive_working,minimum_scba_iterations=2,
        required_consecutive_scba_passes=1,minimum_poisson_iterations=1,
        required_consecutive_poisson_passes=1,stagnation_window=0,stagnation_relative_improvement=0.0,
        diagnostic_quality=DiagnosticQualityPolicy(enabled=true,keldysh_threshold=0.1))
    options=SolverOptions(max_scba=3,max_poisson=3,convergence=policy)
    rows=[SCBAIteration(i,0,0,0,0,0,1,0,0,0,i==1 ? Inf : 0,
        i==1 ? Inf : 0,0,1,1,1,i==1 ? NaN : 0,nothing,nothing) for i in 1:2]
    old=f.solution.scba
    scba=SCBAResult(old.green,old.scattering,old.embedding,old.embedding_plus,old.embedding_minus,
        rows,true,:converged,:strictly_converged,QCLNEGF._solver_restart_contract(options,f.production.algorithms),old.mixer_state)
    outer=[OuterIteration(1,0,0,0,0,0,0,0,0,0,ones(1),ones(25))]
    metrics=Dict(name=>0.0 for name in keys(N._solution_validation_limits(options.tolerances,Dict())))
    approximate && (metrics[:r_K]=0.05)
    report=ConvergenceReport(true,metrics,["manufactured scalar-band control"])
    observations=deepcopy(f.solution.observables)
    observations[:quality]=approximate ? :approximate : :strict
    solution=replace_solution(f.solution;options,scba,outer_history=outer,observables=observations,
        report,converged=!approximate,status=approximate ? :approximate : :converged)
    identity=exact_restart_identity(f.problem,options,f.production,f.execution;
        domain_adaptation=f.policy,kernel_options=nothing)
    state=replace_public_state(f.state;identity,scba,outer_history=outer,status=solution.status,
        progress=SolverProgress(1,1,2,1,1,1),mixer_position=MixerPosition(1,(),0,:none,:none),
        final_evidence=N._final_evidence(solution))
    bind=(;f.binding...,options)
    validate_exact_state(f.problem,state;bind...)
    return (;f...,options,solution,state,binding=bind)
end
# An error receipt must not bless a complete passing scalar-band report. These
# are constructor-only consistency controls; no selected map is evaluated.
@testset "Error revalidation rejects passing reports and preserves failed raw evidence" begin
    for approximate in (false,true)
        f=classified_fixture(;approximate)
        context=AcceptanceContext(backend=:production,algorithms=f.production.algorithms,
            hilbert_roundoff_measured=true,map_kernels_shared_with_solver=true,final_state_alignment=:restored_aligned)
        observations=deepcopy(f.solution.observables)
        observations[:acceptance_measurement_context]=Dict("backend"=>"production",
            "hilbert_roundoff_measured"=>true,"map_kernels_shared_with_solver"=>true,
            "final_state_alignment"=>"restored_aligned")
        observations[:restoration_evidence]=CompletedStateEvidence(original=f.state.final_evidence,
            revalidation=RevalidationReceipt(;status=:error,reason=:selected_map_measurement_failed,
                detail="manufactured error with complete passing report",context,report=f.solution.report),
            independent_oracle=IndependentOracleReceipt(status=:not_measured,
                reason=:not_requested,detail="",context=nothing,report=nothing))
        solution=replace_solution(f.solution;observables=observations)
        before=N._exact_digest((solution,f.state))
        @test caught(()->take_snapshot((;f...,solution);source_mode=:completed_restore)) isa ArgumentError
        @test N._exact_digest((solution,f.state))==before
    end
    # The owning error path may keep converged inner history but publishes a
    # failed current report. Retain every failed IEEE diagnostic unchanged.
    f=classified_fixture()
    payload=reinterpret(Float64,0x7ff8000000000042)
    metrics=copy(f.solution.report.metrics)
    metrics[:r_K]=payload;metrics[:r_Σ]=Inf;metrics[:r_U]=-0.0
    report=ConvergenceReport(false,metrics,["manufactured failed selected-map measurement"])
    context=AcceptanceContext(backend=:production,algorithms=f.production.algorithms,
        hilbert_roundoff_measured=nothing,map_kernels_shared_with_solver=true,final_state_alignment=:restored_aligned)
    observations=deepcopy(f.solution.observables)
    observations[:quality]=:unconverged;observations[:termination_reason]="validation_failed"
    observations[:acceptance_measurement_context]=Dict("backend"=>"production",
        "hilbert_roundoff_measured"=>false,"map_kernels_shared_with_solver"=>true,
        "final_state_alignment"=>"restored_aligned")
    observations[:restoration_evidence]=CompletedStateEvidence(original=f.state.final_evidence,
        revalidation=RevalidationReceipt(;status=:error,reason=:selected_map_measurement_failed,
            detail="manufactured failed raw report",context,report),
        independent_oracle=IndependentOracleReceipt(status=:not_measured,
            reason=:not_requested,detail="",context=nothing,report=nothing))
    solution=replace_solution(f.solution;observables=observations,report,converged=false,status=:validation_failed)
    before=N._exact_digest((solution,f.state))
    captured=take_snapshot((;f...,solution);source_mode=:completed_restore)
    @test captured.revalidation.status===:error
    @test !captured.revalidation.report.passed && !captured.summary.converged && captured.summary.quality===:unconverged
    @test reinterpret(UInt64,captured.revalidation.report.metrics[:r_K])==0x7ff8000000000042
    @test captured.revalidation.report.metrics[:r_Σ]===Inf
    @test signbit(captured.revalidation.report.metrics[:r_U])
    @test f.state.final_evidence.report.passed && !captured.summary.stationary_candidate_accepted
    @test N._exact_digest((solution,f.state))==before
end

@testset "Strict and approximate classifications remain distinct for produced and restored reports" begin
    for originally_approximate in (false,true)
        f=classified_fixture(approximate=originally_approximate)
        produced=take_snapshot(f)
        @test produced.summary.converged==!originally_approximate
        @test produced.summary.quality== (originally_approximate ? :approximate : :strict)
        @test produced.summary.scientific_accepted===false
        target=classified_fixture(approximate=!originally_approximate)
        context=AcceptanceContext(backend=:production,algorithms=f.production.algorithms,
            hilbert_roundoff_measured=true,map_kernels_shared_with_solver=true,final_state_alignment=:restored_aligned)
        report=target.solution.report
        observation=deepcopy(f.solution.observables)
        observation[:quality]=target.solution.observables[:quality]
        observation[:acceptance_measurement_context]=Dict("backend"=>"production",
            "hilbert_roundoff_measured"=>true,"map_kernels_shared_with_solver"=>true,
            "final_state_alignment"=>"restored_aligned")
        observation[:restoration_evidence]=CompletedStateEvidence(original=f.state.final_evidence,
            revalidation=RevalidationReceipt(;status=:measured,reason=nothing,detail="constructor control",
                context,report),independent_oracle=IndependentOracleReceipt(status=:not_measured,
                reason=:not_requested,detail="",context=nothing,report=nothing))
        solution=replace_solution(f.solution;report,observables=observation,
            converged=target.solution.converged,status=target.solution.status)
        restored=take_snapshot((;f...,solution);source_mode=:completed_restore)
        @test restored.summary.converged==originally_approximate
        @test restored.summary.quality==(originally_approximate ? :strict : :approximate)
        @test restored.independent_oracle.status===:not_measured
        @test restored.summary.scientific_accepted===false
        @test_throws ArgumentError CoreRecord(Dict("snapshot"=>restored))
        wrong=replace_solution(solution;converged=!solution.converged)
        @test caught(()->take_snapshot((;f...,solution=wrong);source_mode=:completed_restore)) isa ArgumentError
    end
end

@testset "Legacy effective views use the first-k upper Hermitian triangle, not unrelated samples" begin
    f=final_fixture(basis_states=2)
    f.problem.grids.κ[2]=NaN
    f=reseal(f)
    capture=take_snapshot(f)
    @test capture.views.effective_levels_eV!==nothing
    @test capture.views.effective_levels_eV!==nothing && all(isfinite,capture.views.effective_levels_eV)
    @test isnan(capture.context.momentum_scaled[2])
    f=final_fixture(basis_states=2)
    f.problem.basis.H₀[2,1]=complex(NaN,0.0)
    f=reseal(f)
    capture=take_snapshot(f)
    @test capture.views.effective_transform!==nothing
    @test capture.views.effective_levels_eV!==nothing && all(isfinite,capture.views.effective_levels_eV)
    f=final_fixture()
    f.solution.Uᴴ[1]=NaN # parabolic owner returns raw projected NaN values
    f=reseal(f)
    capture=take_snapshot(f)
    @test capture.views.localized_levels_eV!==nothing
    @test capture.views.localized_levels_eV!==nothing && any(isnan,capture.views.localized_levels_eV)
    @test capture.views.effective_transform===nothing
end

@testset "Nonempty Anderson payload participates in exact observation binding" begin
    f=final_fixture()
    production=QCLNEGF.with_production_options(f.production;
        algorithms=AlgorithmOptions(mixing=:anderson,anderson_history_depth=2))
    execution=public_state_execution(production)
    old=f.solution.scba
    names=fieldnames(SCBAIteration)
    rows=[SCBAIteration((name===:ν ? i : getfield(old.history[1],name) for name in names)...) for i in 1:3]
    states=[[deepcopy(old.embedding_plus),deepcopy(old.embedding_minus)] for _ in 1:2]
    residuals=deepcopy(states)
    states[1][1].Σᴿ[1]=0.25+0.5im
    mixer=SCBAMixerState(:anderson,states,residuals)
    scba=SCBAResult(old.green,old.scattering,old.embedding,old.embedding_plus,old.embedding_minus,
        rows,false,:invalid_candidate,:invalid,QCLNEGF._solver_restart_contract(f.options,production.algorithms),mixer)
    solution=replace_solution(f.solution;scba)
    identity=exact_restart_identity(f.problem,f.options,production,execution;
        domain_adaptation=f.policy,kernel_options=nothing)
    state=replace_public_state(f.state;identity,scba,progress=SolverProgress(0,1,3,0,0,0),
        mixer_position=MixerPosition(2,(1,2),0,:none,:none),
        measurement_context=AcceptanceContext(backend=:production,algorithms=production.algorithms,
            hilbert_roundoff_measured=true,map_kernels_shared_with_solver=true,final_state_alignment=:producer_aligned),
        final_evidence=N._final_evidence(solution))
    binding=(;f.binding...,production_options=production,execution)
    @test validate_exact_state(f.problem,state;binding...)===nothing
    packet=(;f...,solution,state,binding)
    @test take_snapshot(packet).terminal_integrity==state.integrity
    wrong=deepcopy(solution);wrong.scba.mixer_state.states[1][1].Σᴿ[1]+=im
    @test caught(()->take_snapshot((;packet...,solution=wrong))) isa ArgumentError
    wrong=deepcopy(solution);empty!(wrong.scba.mixer_state.residuals)
    @test caught(()->take_snapshot((;packet...,solution=wrong))) isa ArgumentError
end
end
