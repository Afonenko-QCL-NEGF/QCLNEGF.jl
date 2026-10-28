function public_state_problem(;scattering_options=ScatteringOptions(
    LO=false,acoustic=false,impurity=false,IFR=false,alloy=false),
    physical_models=PhysicalModelOptions())
    base=tutorial_numerics()
    changes=(N_z=25,N_b=1,P_basis=1,N_E=9,N_k=2,N_qz=3,N_φ=4)
    numerical=NumericalParameters(;
        (name=>(hasproperty(changes,name) ? getproperty(changes,name) : getfield(base,name))
         for name in fieldnames(NumericalParameters))...)
    # Tiny synthetic domain: explicitly reduce ionized target instead of bypassing Pauli capacity.
    return build_problem(physical=reference_parameters(f_ion=1/450),
        numerical=numerical,scattering=scattering_options,physical_models=physical_models)
end
function public_state_execution(production)
    return ExactExecutionIdentity(core_code="p02-synthetic-source",julia=string(VERSION),
        blas=string(BLAS.get_config()),fftw=string(Base.pkgversion(QCLNEGF.FFTW)),
        scalar_type="Float64/ComplexF64",index_type=string(Int),reduction="static_chunks_v1",
        worker_count=production.worker_count,blas_threads=BLAS.get_num_threads(),
        energy_chunk=production.energy_chunk,hilbert_columns=production.hilbert_columns,
        residual_chunk=production.residual_chunk,parallel_backend=production.parallel_backend)
end
function public_state_fixture(;scattering_options=ScatteringOptions(
    LO=false,acoustic=false,impurity=false,IFR=false,alloy=false),
    physical_models=PhysicalModelOptions())
    problem=public_state_problem(;scattering_options,physical_models)
    options=SolverOptions(max_scba=3,max_poisson=3)
    production=ProductionOptions(worker_count=1,checkpoint_every_scba=1,checkpoint_every_outer=1,
        progress_every_scba=0,progress_every_outer=0)
    execution=public_state_execution(production)
    policy=DomainAdaptationPolicy()
    identity=exact_restart_identity(problem,options,production,execution;
        domain_adaptation=policy,kernel_options=nothing)
    shape=(9,2,1,1)
    green=GreenState(fill(1.0-1im,shape),fill(0.5im,shape),fill(-1.5im,shape),
        fill(2.0+0im,shape),ones(9,2),ones(9,2))
    family()=SelfEnergyFamily(zeros(ComplexF64,shape),zeros(ComplexF64,shape),zeros(ComplexF64,shape))
    history=[SCBAIteration(1,0.2,0.3,0.4,0.5,0.0,1.0,0.0,0.0,0.0,
        Inf,Inf,3.0,1.0,1.0,1.0,NaN,nothing,nothing)]
    scattering=Dict{Symbol,SelfEnergyFamily}(name=>family() for name in problem.kernels.enabled)
    scba=SCBAResult(green,scattering,family(),family(),family(),history,
        false,:running,:unresolved,
        QCLNEGF._solver_restart_contract(options,production.algorithms),SCBAMixerState())
    continuation=SCBAContinuation(decision=SCBAStopDecision(false,false,:running,:unresolved,"pending"),
        candidate=FixedPointAuditCandidate(Dict{Symbol,SelfEnergyFamily}(
            name=>family() for name in problem.kernels.enabled),family(),family(),family(),0.0),
        next_previous_current=3.0,next_previous_populations=
            QCLNEGF.QCLNumerics._production_state_populations(problem,green.Gˡ),seed_mu_scaled=nothing,
        seed_number_ratio=nothing,witness_retention=(0,0,0))
    state=ExactSolverState(identity=identity,boundary=:scba_accepted,published_U=zeros(25),
        published_n=ones(25),scba=scba,mixer_position=QCLNEGF.QCLNumerics._initial_mixer_position(),
        family_order=Tuple((problem.kernels.enabled...,:embedding_plus,:embedding_minus)),outer_history=OuterIteration[],
        progress=SolverProgress(completed_outer=0,current_outer=1,accepted_scba=1,
            outer_strict_streak=0,inner_strict_streak=0,diagnostic_streak=0),
        continuation=continuation,status=:running_scba,
        warnings=(CoreRecord(Dict("code"=>"SYNTHETIC","unmeasured"=>NaN)),),policy_evaluations=(),adaptation=nothing,
        measurement_context=AcceptanceContext(backend=:production,algorithms=production.algorithms,
            hilbert_roundoff_measured=true,map_kernels_shared_with_solver=true,
            final_state_alignment=:producer_aligned),final_evidence=nothing)
    return (;problem,options,production,execution,policy,state)
end

function public_state_failed_scba(old)
    history=deepcopy(old.history)
    row=last(history)
    names=fieldnames(SCBAIteration)
    values=NamedTuple{names}(Tuple(getfield(row,name) for name in names))
    changed=merge(values,(;r_K=Inf,r_Σ=Inf))
    history[end]=SCBAIteration((getfield(changed,name) for name in names)...)
    return SCBAResult(old.green,old.scattering,old.embedding,old.embedding_plus,old.embedding_minus,
        history,false,:invalid_candidate,:invalid,old.restart_contract,old.mixer_state)
end

function replace_public_state(state;changes...)
    names=fieldnames(ExactSolverState)[1:end-1]
    values=NamedTuple{names}(Tuple(getfield(state,name) for name in names))
    return ExactSolverState(;merge(values,(;changes...))...)
end
