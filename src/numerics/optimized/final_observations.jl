# Owning final postprocessing only. No file/transport code or new physical gates.
import ..QCLPhysics: spatial_energy_density
_final_add(a::Int,b::Int)=Base.Checked.checked_add(a,b)
_final_mul(a::Int,b::Int)=Base.Checked.checked_mul(a,b)
_final_sum(values)=foldl(_final_add,values;init=0)
_final_require(ok,path)=ok || throw(ArgumentError("final observation binding: $path"))


function _final_check_context(;dimensions,domain_revision,discretization_id,basis_id,
    x_scaled,wX,energy_scaled,wE,momentum_scaled,wk,z_nm,energy_eV,k_per_nm,
    state_index,trusted_energy,kwargs...)
    z,e,k,b=dimensions
    all(>(0),dimensions) || throw(DimensionMismatch("final context has empty axes"))
    for (a,n) in ((x_scaled,z),(wX,z),(z_nm,z),(energy_scaled,e),(wE,e),(energy_eV,e),
                  (momentum_scaled,k),(wk,k),(k_per_nm,k),(state_index,b),(trusted_energy,e))
        length(a)==n || throw(DimensionMismatch("final context coordinate/weight shape differs"))
    end
    _final_require(domain_revision>=0,"domain revision is negative")
    _final_require(all(i->state_index[i]==i-1,eachindex(state_index)),"state coordinates differ")
    _final_require(all(id->occursin(r"^[0-9a-f]{64}\z",id),(discretization_id,basis_id)),"domain IDs differ")
    nothing
end
function _final_check_views(;potential_structure_eV,potential_external_eV,potential_hartree_eV,
    potential_total_eV,density_per_m3,spatial_energy_density_per_eV_m3,spectral_trace_per_eV,
    occupied_trace_per_eV,localized_levels_eV,effective_levels_eV,localized_wavefunctions_per_sqrt_nm,
    effective_wavefunctions_per_sqrt_nm,effective_transform,sheet_density_matrix_per_m2,
    effective_density_matrix_per_m2,kwargs...)
    z=length(potential_structure_eV);e,k=size(spectral_trace_per_eV);b=size(sheet_density_matrix_per_m2,1)
    for array in (potential_external_eV,potential_hartree_eV,potential_total_eV,density_per_m3)
        length(array)==z || throw(DimensionMismatch("final view spatial axis differs"))
    end
    size(occupied_trace_per_eV)==(e,k) || throw(DimensionMismatch("final trace shape differs"))
    for (array,shape) in ((spatial_energy_density_per_eV_m3,(e,z)),
        (localized_wavefunctions_per_sqrt_nm,(z,b)),(effective_wavefunctions_per_sqrt_nm,(z,b)),
        (localized_levels_eV,(b,)),(effective_levels_eV,(b,)),(effective_transform,(b,b)),
        (sheet_density_matrix_per_m2,(b,b)),(effective_density_matrix_per_m2,(b,b)))
        array===nothing || size(array)==shape || throw(DimensionMismatch("final view basis/energy shape differs"))
    end
    group=(effective_levels_eV,effective_wavefunctions_per_sqrt_nm,effective_transform,effective_density_matrix_per_m2)
    _final_require(all(isnothing,group) || all(!isnothing,group),"incomplete effective view group")
    nothing
end
function _final_check_snapshot(;source_mode,context,views,positivity,revalidation,independent_oracle,issues,kwargs...)
    _final_require(source_mode in (:produced,:completed_restore),"unknown source mode")
    _final_require(source_mode===:produced ? revalidation===nothing && independent_oracle===nothing :
        revalidation!==nothing && independent_oracle!==nothing,"receipt mode differs")
    z,e,k,b=context.dimensions
    size(views.spectral_trace_per_eV)==(e,k) && length(views.density_per_m3)==z &&
        size(views.sheet_density_matrix_per_m2)==(b,b) || throw(DimensionMismatch("context/view axes differ"))
    targets=Tuple(i.target for i in issues)
    _final_require(length(unique(targets))==length(targets),"duplicate observation issues")
    for name in fieldnames(StoredFinalObservations)
        if getfield(views,name)===nothing
            _final_require(name in targets,"unavailable $name has no issue")
        end
    end
    _final_require(positivity!==nothing || :positivity in targets,"unavailable positivity has no issue")
    if positivity!==nothing && positivity.integrals!==nothing
        integral=positivity.integrals
        size(integral.spectral_integral)==(k,b,b) && size(integral.spectral_integral_eigenvalues)==(k,b) &&
            size(integral.spectral_trace_cumulative)==(e,k) && size(integral.negative_energy_integral_scaled)==(6,k) ||
            throw(DimensionMismatch("final integral/context axes differ"))
    end
    nothing
end

function _final_dimensions(solution,terminal)
    p=solution.problem; n=p.numerical; g=p.grids; b=p.basis
    dims=(n.N_z,n.N_E,n.N_k,n.N_b)
    all(>(0),dims) || throw(DimensionMismatch("final observation axes must be nonempty"))
    z,e,k,nb=dims
    for (value,count) in ((g.x,z),(g.wˣ,z),(g.ε,e),(g.wᴱ,e),(g.trusted_energy,e),
                         (g.κ,k),(g.wᵏ,k),(solution.Uᴴ,z),(solution.n,z),(p.profiles.Eᶜ,z))
        length(value)==count || throw(DimensionMismatch("final observation coordinate/profile axis differs"))
    end
    size(b.χ)==size(b.Φ)==(z,nb) || throw(DimensionMismatch("final observation basis shape differs"))
    for value in (b.H₀,b.Z,b.M⁻¹)
        size(value)==(nb,nb) || throw(DimensionMismatch("final observation projected operator differs"))
    end
    _check_diagnostic_state_shapes(solution)
    return dims
end

function _final_numeric_sizes(dimensions::NTuple{4,Int})
    z,e,k,b=dimensions
    zb,bb,ek,ez,kbb,kb=(_final_mul(z,b),_final_mul(b,b),_final_mul(e,k),
        _final_mul(e,z),_final_mul(k,_final_mul(b,b)),_final_mul(k,b))
    c=_final_add(_final_mul(8,_final_sum((_final_mul(3,z),_final_mul(3,e),_final_mul(3,k),b))),max(e,_final_mul(sizeof(UInt64),cld(e,64))))
    v=_final_sum((_final_mul(8,_final_sum((_final_mul(5,z),ez,_final_mul(2,ek),_final_mul(2,b)))),
        _final_mul(16,_final_sum((_final_mul(2,zb),_final_mul(3,bb)))),16))
    d=_final_add(_final_mul(16,kbb),_final_mul(8,_final_sum((kb,ek,_final_mul(6,k)))))
    a=_final_mul(8,_final_sum((z,e,k)))
    w=_final_sum((_final_mul(16,_final_sum((kbb,_final_mul(8,zb),_final_mul(64,bb)))),
        _final_mul(8,_final_sum((_final_mul(4,z),_final_mul(4,e),_final_mul(4,k),_final_mul(64,b)))),65536))
    return (;context=c,views=v,positivity=d,assessment=a,workspace=w)
end

# Julia 1.13 Complex Hermitian eigen/eigvals default to LAPACK ZHEEVR.
# Workspace-query contract: LAPACK 3.12.0 SRC/zheevr.f returns at LQUERY,
# before W/Z/ISUPPZ/A access. Only WORK(1), RWORK(1), IWORK(1), INFO are written.
# https://github.com/Reference-LAPACK/lapack/blob/v3.12.0/SRC/zheevr.f
const _FINAL_ZHEEVR=LinearAlgebra.BLAS.@blasfunc(zheevr_)
const _FINAL_BLAS_INT=LinearAlgebra.BlasInt
function _final_query_count(value)
    isfinite(value) && value>=1 && isinteger(value) ||
        throw(ArgumentError("invalid LAPACK workspace query count"))
    return Int(value) # checked representability, including f64 values above typemax(Int)
end
function _final_eigen_workspace_bytes(b::Int)
    b>0 || throw(ArgumentError("LAPACK workspace dimension must be positive"))
    # These driver dimensions/products must fit the actual Fortran integer ABI,
    # before a foreign call or a potentially overflowing LAPACK integer product.
    for value in (b,_final_mul(b,b),_final_mul(24,b),_final_mul(10,b))
        value<=typemax(_FINAL_BLAS_INT) || throw(OverflowError("LAPACK dimension is not representable"))
    end
    _final_mul(16,_final_mul(b,b))
    n=_FINAL_BLAS_INT(b)
    maximum=0
    for job in (UInt8('V'),UInt8('N'))
        a=Ref{ComplexF64}(complex(NaN,NaN));w=Ref{Float64}(91.0);z=Ref{ComplexF64}(92.0+3im)
        support=Ref{_FINAL_BLAS_INT}(93);m=Ref{_FINAL_BLAS_INT}(94)
        work=Ref{ComplexF64}(0);rwork=Ref{Float64}(0);iwork=Ref{_FINAL_BLAS_INT}(0)
        info=Ref{_FINAL_BLAS_INT}(0)
        ccall((_FINAL_ZHEEVR,LinearAlgebra.LAPACK.libblastrampoline),Cvoid,
            (Ref{UInt8},Ref{UInt8},Ref{UInt8},Ref{_FINAL_BLAS_INT},
             Ref{ComplexF64},Ref{_FINAL_BLAS_INT},Ref{Float64},Ref{Float64},
             Ref{_FINAL_BLAS_INT},Ref{_FINAL_BLAS_INT},Ref{Float64},Ref{_FINAL_BLAS_INT},
             Ref{Float64},Ref{ComplexF64},Ref{_FINAL_BLAS_INT},Ref{_FINAL_BLAS_INT},
             Ref{ComplexF64},Ref{_FINAL_BLAS_INT},Ref{Float64},Ref{_FINAL_BLAS_INT},
             Ref{_FINAL_BLAS_INT},Ref{_FINAL_BLAS_INT},Ref{_FINAL_BLAS_INT},Clong,Clong,Clong),
            job,UInt8('A'),UInt8('U'),n,a,n,0.0,0.0,
            _FINAL_BLAS_INT(0),_FINAL_BLAS_INT(0),-1.0,m,w,z,
            job==UInt8('V') ? n : _FINAL_BLAS_INT(1),support,
            work,_FINAL_BLAS_INT(-1),rwork,_FINAL_BLAS_INT(-1),
            iwork,_FINAL_BLAS_INT(-1),info,Clong(1),Clong(1),Clong(1))
        info[]==0 || throw(ArgumentError("LAPACK workspace query failed: INFO=$(info[])"))
        # Sentinel outputs are not part of a workspace query. Never retry as an eigensolve.
        _final_require(w[]==91.0 && z[]==92.0+3im && support[]==93 && m[]==94 &&
            isnan(real(a[])) && isnan(imag(a[])),"LAPACK query accessed nonworkspace outputs")
        lwork=_final_query_count(real(work[]));lrwork=_final_query_count(rwork[])
        liwork=_final_query_count(iwork[])
        # WORK/RWORK/IWORK, returned W, support indices. Eigenvector output and
        # input matrix copies are already in V / numeric baseline W respectively.
        bytes=_final_sum((_final_mul(16,lwork),_final_mul(8,lrwork),
            _final_mul(sizeof(_FINAL_BLAS_INT),liwork),_final_mul(8,b),
            _final_mul(2sizeof(_FINAL_BLAS_INT),b)))
        maximum=max(maximum,bytes)
    end
    return maximum
end

"""
    final_observation_memory_estimate(solution, terminal; maximum_metadata_bytes)

Check shapes and checked-Int payload sizes before dominant allocations. Counts
are conservative array/metadata allowances, not total process RSS. Existing
canonical validation buffers are explicit. A scalar LAPACK workspace query is
performed; no matrix allocation, eigensolve, getter, hash or map is executed.
"""
function final_observation_memory_estimate(solution::NEGFSolution,terminal::ExactSolverState;
    maximum_metadata_bytes::Int)
    maximum_metadata_bytes>0 || throw(ArgumentError("maximum_metadata_bytes must be positive"))
    dimensions=_final_dimensions(solution,terminal)
    c,v,d,a,w=values(_final_numeric_sizes(dimensions))
    # Only data which the output actually detaches: no untyped Green/map cache.
    audit=get(solution.observables,:final_fixed_point_audit,nothing)
    receipt=get(solution.observables,:restoration_evidence,nothing)
    metadata=Base.summarysize((solution.options,solution.report,audit,receipt,
        terminal.final_evidence,terminal.warnings,solution.problem.models,solution.problem.scattering))
    m=_final_add(65536,_final_mul(8,metadata))
    m<=maximum_metadata_bytes || throw(ArgumentError("final observation metadata exceeds maximum_metadata_bytes"))
    w=_final_add(w,_final_eigen_workspace_bytes(dimensions[4]))
    baseline=_final_add(_final_mul(4,Base.summarysize((solution.problem,terminal))),32768)
    residency=Base.summarysize((solution,terminal))
    return FinalObservationMemoryEstimate(dimensions,residency,c,v,d,m,w,baseline,
        _final_sum((c,v,d,a,m,w,baseline)))
end

# IEEE payload comparison, not isequal (which identifies distinct NaNs/zero signs).
# Only the finite owning public payload closure is traversed; no serialization or copy.
function _final_same_bits(a,b)
    typeof(a)===typeof(b) || return false
    if a isa Float64
        return reinterpret(UInt64,a)==reinterpret(UInt64,b)
    elseif a isa Float32
        return reinterpret(UInt32,a)==reinterpret(UInt32,b)
    elseif a isa Complex
        return _final_same_bits(real(a),real(b)) && _final_same_bits(imag(a),imag(b))
    elseif a isa Union{Nothing,Bool,Integer,Symbol,AbstractString}
        return a==b
    elseif a isa AbstractArray
        axes(a)==axes(b) || return false
        return all(i->_final_same_bits(a[i],b[i]),eachindex(a,b))
    elseif a isa AbstractDict
        length(a)==length(b) || return false
        return all(key->haskey(b,key) && _final_same_bits(a[key],b[key]),keys(a))
    elseif a isa Tuple
        return length(a)==length(b) && all(i->_final_same_bits(a[i],b[i]),eachindex(a))
    end
    isstructtype(typeof(a)) || return false
    return all(name->_final_same_bits(getfield(a,name),getfield(b,name)),fieldnames(typeof(a)))
end
function _final_core_matches(saved::CoreScalar,value)
    _final_same_bits(saved.value,value)
end
function _final_core_matches(saved::CoreArray,value)
    value isa AbstractArray && size(value)==saved.shape && eltype(value)===eltype(saved.values) || return false
    all(i->_final_same_bits(saved.values[i],value[i]),eachindex(saved.values))
end
function _final_core_matches(saved::CoreList,value)
    value isa Union{AbstractVector,Tuple} && length(value)==length(saved.values) || return false
    all(i->_final_core_matches(saved.values[i],value[i]),eachindex(saved.values))
end
function _final_core_matches(saved::CoreRecord,value)
    value isa AbstractDict && length(value)==length(saved.entries) || return false
    all(saved.entries) do pair
        key=first(pair);actual=haskey(value,key) ? key : Symbol(key)
        haskey(value,actual) && _final_core_matches(last(pair),value[actual])
    end
end

# Reuse the owning scalar band evaluator, not a fresh map or a copied threshold.
# A typed receipt/report alone cannot manufacture strict acceptance against the
# preserved inner/outer history. Original and revalidated reports stay separate.
function _final_restored_classification(solution,state)
    strict_report=_validate_solution_metrics(solution,solution.report.metrics)
    outer=state.progress.outer_strict_streak>=solution.options.convergence.required_consecutive_poisson_passes
    strict=outer && solution.scba.converged && strict_report.passed
    approximate=false
    if !strict && outer && scba_accepted(solution.scba) &&
        solution.options.convergence.mode===:adaptive_working &&
        solution.options.convergence.diagnostic_quality.enabled
        band=NEGFSolution(solution.problem,_approximate_options(solution.options),solution.Uᴴ,
            solution.n,solution.scba,solution.outer_history,solution.observables,solution.report,
            false,solution.status)
        approximate=_validate_solution_metrics(band,solution.report.metrics;accept_approximate_scba=true).passed
    end
    _final_require(solution.report.passed==(strict_report.passed || approximate),
        "restored report contradicts owning scalar band gates")
    expected=strict ? :converged : approximate ? :approximate :
        solution.report.passed ? state.status : :validation_failed
    _final_require(solution.converged==strict && solution.status===expected,
        "restored classification contradicts history/report")
    nothing
end

function _final_bind(solution,state,options,production,source_mode)
    _final_require(state.boundary===:final && state.final_evidence!==nothing,"global final is required")
    _final_require(_final_same_bits(solution.options,options),"options differ")
    for (name,a,b) in ((:U,solution.Uᴴ,state.published_U),(:n,solution.n,state.published_n),
        (:scba,solution.scba,state.scba),(:outer_history,solution.outer_history,state.outer_history))
        _final_require(_final_same_bits(a,b),"$name differs from terminal")
    end
    warnings=get(solution.observables,:warnings,nothing)
    _final_require(warnings isa AbstractVector && length(warnings)==length(state.warnings) &&
        all(i->_final_core_matches(state.warnings[i],warnings[i]),eachindex(state.warnings)),"warnings differ")
    if source_mode===:produced
        _final_require(!haskey(solution.observables,:restoration_evidence),"unexpected restoration evidence")
        ev=state.final_evidence
        _final_require(_final_same_bits(solution.report,ev.report) &&
            solution.converged==ev.original_converged && solution.status==ev.original_status,"original classification differs")
        actual=Dict{Symbol,Any}(key=>solution.observables[key] for key in _EXACT_FINAL_OBSERVATIONS
            if haskey(solution.observables,key))
        _final_require(_final_core_matches(ev.observations,actual),"original observations differ")
        _final_require(_final_same_bits(_acceptance_context(production;observations=solution.observables),
            state.measurement_context),"measurement context differs")
        return nothing,nothing
    end
    ev=get(solution.observables,:restoration_evidence,nothing)
    _final_require(ev isa CompletedStateEvidence,"completed restoration receipt is missing")
    _final_require(_final_same_bits(ev.original,state.final_evidence),"restoration original differs")
    _final_require(_final_same_bits(solution.report,ev.revalidation.report),"revalidated report differs")
    ctx=ev.revalidation.context
    _final_require(_final_same_bits(ctx.algorithms,production.algorithms) &&
        ctx.final_state_alignment===:restored_aligned,"restoration context differs")
    raw=get(solution.observables,:acceptance_measurement_context,nothing)
    _final_require(raw isa AbstractDict,"restoration measurement context is missing")
    for (key,expected) in (("backend",String(ctx.backend)),
        ("map_kernels_shared_with_solver",ctx.map_kernels_shared_with_solver),
        ("final_state_alignment",String(ctx.final_state_alignment)))
        _final_require(haskey(raw,key) && _final_same_bits(raw[key],expected),"restoration $key differs")
    end
    measured=get(raw,"hilbert_roundoff_measured",nothing)
    _final_require(_final_same_bits(measured,ctx.hilbert_roundoff_measured) ||
        (ev.revalidation.status===:error && ctx.hilbert_roundoff_measured===nothing && measured===false),
        "restoration Hilbert measurement differs")
    # A failed measurement cannot carry a passing report. The existing scalar
    # band/history check below also binds the current converged/status fields.
    _final_require(ev.revalidation.status!==:error || !solution.report.passed,
        "error revalidation has a passing report")
    _final_restored_classification(solution,state)
    if ev.revalidation.status===:measured
        _final_require(ctx.backend===:production && ctx.map_kernels_shared_with_solver===true,
            "measured revalidation lacks selected-map context")
    end
    independent=ev.independent_oracle
    if independent.status===:measured
        _final_require(independent.report!==nothing && independent.context!==nothing &&
            independent.context.backend===:reference && independent.context.map_kernels_shared_with_solver===false,
            "measured independent oracle lacks independent context/report")
    end
    return deepcopy(ev.revalidation),deepcopy(ev.independent_oracle)
end

function _final_context(solution,state,production)
    p=solution.problem;g=p.grids;b=p.basis;s=p.scales;physical=p.physical
    field=Float64(ustrip(u"V/m",physical.F_bias));period=_metres(period_length(physical))
    # Adaptive final validation already verifies these cached receipt IDs against
    # this exact domain. Fixed-domain states need one explicit representation digest.
    domain_ids=if state.adaptation===nothing
        (_representation_digest(g.x,g.wˣ,g.ε,g.wᴱ,g.κ,g.wᵏ,g.qᶻ,g.φ),
         _representation_digest(b.Φ,b.H₀,b.M⁻¹))
    else
        evidence=last(state.adaptation.revisions).evidence
        (evidence.discretization_id,evidence.basis_id)
    end
    context=FinalObservationContext((p.numerical.N_z,p.numerical.N_E,p.numerical.N_k,p.numerical.N_b),
        state.adaptation===nothing ? 0 : state.adaptation.expansions,
        domain_ids...,copy(g.x),copy(g.wˣ),copy(g.ε),copy(g.wᴱ),
        copy(g.κ),copy(g.wᵏ),g.x .* s.L₀_m .* 1e9,g.ε .* s.E₀_eV .+ _electronvolts(physical.E_ref),
        g.κ ./ s.L₀_m .* 1e-9,collect(Int64(0):Int64(p.numerical.N_b-1)),copy(g.trusted_energy),
        FinalObservationScales(s.E₀_eV,s.L₀_m,Float64(ustrip(u"A/m^2",s.J₀))),
        FinalObservationInputs(field,field*period,Float64(ustrip(u"K",physical.Tᴸ)),
            Float64(ustrip(u"K",physical.Tᴸᴼ)),_electronvolts(physical.E_ref),
            Float64(ustrip(u"m^-2",physical.N_dop²ᴰ)),physical.f_ion,_metres(physical.z₀),period,physical.g_s),
        p.scattering,Tuple(name for name in p.kernels.enabled if name in _SCATTERING_MECHANISM_ORDER),
        Tuple(p.kernels.enabled),p.models,production.algorithms,p.energy_shift_discretization,b.localization,
        :stationary_field_periodic_embedding,:orthonormal_finite_projected,:scaled_E0_plus_E_ref,
        :k0_ascending_spatial_first_max_real_positive)
    return context
end
function _final_issue!(issues,target,reason,detail)
    push!(issues,FinalObservationIssue(target,:error,reason,detail))
    nothing
end
_final_expected_error(error)=error isa Union{LinearAlgebra.LAPACKException,LinearAlgebra.SingularException,DomainError}

function _final_views(solution,context,issues)
    p=solution.problem;g=p.grids;s=p.scales;green=solution.scba.green
    structure=p.profiles.Eᶜ .* s.E₀_eV .+ context.inputs.E_ref_eV
    external=-context.inputs.F_bias_V_per_m .* (g.x .* s.L₀_m .- context.inputs.z0_m)
    hartree=solution.Uᴴ .* s.E₀_eV
    density=solution.n ./ s.L₀_m^3
    spectral=[real(sum(green.A[e,k,a,a] for a in axes(green.A,3)))/s.E₀_eV
        for e in eachindex(g.ε),k in eachindex(g.κ)]
    occupied=[real(sum(-im*green.Gˡ[e,k,a,a] for a in axes(green.Gˡ,3)))/s.E₀_eV
        for e in eachindex(g.ε),k in eachindex(g.κ)]
    # This getter explicitly rejects nonfinite GL; retain the raw state/traces instead.
    energy_density=if all(isfinite,green.Gˡ)
        spatial_energy_density(p,green.Gˡ).n_per_eV_m3
    else
        _final_issue!(issues,:spatial_energy_density_per_eV_m3,:nonfinite_green,
            "spatial_energy_density requires finite lesser Green values")
    end
    sheet=ustrip.(u"m^-2",sheet_density_matrix(p,green.Gˡ))
    localized_levels=effective_levels=effective_wavefunctions=transform=effective_sheet=nothing
    projection_reason=:nonfinite_hartree
    projection_detail="Hartree field is nonfinite; projected levels were not evaluated"
    h0=if p.models.dispersion===:parabolic || all(isfinite,solution.Uᴴ)
        try
            h=project_hamiltonians(p,solution.Uᴴ)
            Matrix(view(h,1,:,:))
        catch error
            _final_expected_error(error) || rethrow()
            projection_reason=:projection_failed
            projection_detail=sprint(showerror,error)
            nothing
        end
    else
        nothing
    end
    if h0!==nothing
        # Retain obtained IEEE diagonal samples even when the eigensystem is unavailable.
        localized_levels=real.(diag(h0)) .* s.E₀_eV .+ context.inputs.E_ref_eV
        projection_reason=:nonfinite_hamiltonian
        projection_detail="first-k Hermitian Hamiltonian is nonfinite"
    end
    # Preserve the old upper-Hermitian, first-k domain exactly. Other momentum
    # blocks or the unused lower triangle must not introduce a new rejection gate.
    if h0!==nothing && all(isfinite,Hermitian(h0))
        try
            decomposition=eigen(Hermitian(h0))
            transform=decomposition.vectors
            effective_wavefunctions=p.basis.χ*transform
            for state in axes(effective_wavefunctions,2)
                pivot=argmax(abs.(view(effective_wavefunctions,:,state)))
                amplitude=effective_wavefunctions[pivot,state]
                if abs(amplitude)>0
                    phase=conj(amplitude)/abs(amplitude)
                    effective_wavefunctions[:,state].*=phase
                    transform[:,state].*=phase
                end
            end
            effective_levels=decomposition.values .* s.E₀_eV .+ context.inputs.E_ref_eV
            effective_sheet=transform'*sheet*transform
            effective_wavefunctions ./= sqrt(s.L₀_m*1e9)
        catch error
            _final_expected_error(error) || rethrow()
            effective_levels=effective_wavefunctions=transform=effective_sheet=nothing
            for target in (:effective_levels_eV,:effective_wavefunctions_per_sqrt_nm,
                           :effective_transform,:effective_density_matrix_per_m2)
                _final_issue!(issues,target,:eigensystem_failed,sprint(showerror,error))
            end
        end
    else
        for target in (:effective_levels_eV,:effective_wavefunctions_per_sqrt_nm,
                       :effective_transform,:effective_density_matrix_per_m2)
            _final_issue!(issues,target,projection_reason,projection_detail)
        end
        h0===nothing && _final_issue!(issues,:localized_levels_eV,projection_reason,projection_detail)
    end
    current=Float64(ustrip(u"A/m^2",CODATA.e*boundary_flux(p,solution.scba.embedding_plus,green)))
    return StoredFinalObservations(structure,external,hartree,structure .+ external .+ hartree,
        density,energy_density,spectral,occupied,localized_levels,effective_levels,
        p.basis.χ ./ sqrt(s.L₀_m*1e9),effective_wavefunctions,transform,sheet,effective_sheet,
        current,real(tr(sheet)))
end

"""
    final_observation_snapshot(solution, terminal; options, production_options,
        execution, domain_adaptation, kernel_options, source_mode,
        maximum_bytes, maximum_metadata_bytes)

Capture detached actual context, legacy stored views and one integrated PSD scan.
The supplied global-final snapshot is validated once and compared bitwise to the
solution before physical postprocessing. Explicit :produced/:completed_restore
policies separate original evidence from existing restoration receipts. No fresh
map, solver, hidden restore, threshold repair or file I/O is performed.
"""
function final_observation_snapshot(solution::NEGFSolution,terminal::ExactSolverState;
    options::SolverOptions,production_options::ProductionOptions,
    execution::ExactExecutionIdentity,domain_adaptation::DomainAdaptationPolicy,
    kernel_options::Union{Nothing,ProductionKernelOptions},source_mode::Symbol,
    maximum_bytes::Int,maximum_metadata_bytes::Int)
    source_mode in (:produced,:completed_restore) || throw(ArgumentError("unknown final observation source_mode"))
    maximum_bytes>0 || throw(ArgumentError("maximum_bytes must be positive"))
    _final_require(terminal.boundary===:final && terminal.final_evidence!==nothing,"global final is required")
    estimate=final_observation_memory_estimate(solution,terminal;maximum_metadata_bytes)
    estimate.peak_additional_bytes<=maximum_bytes || throw(ArgumentError("final observations exceed maximum_bytes"))
    validate_exact_state(solution.problem,terminal;options,production_options,execution,
        domain_adaptation,kernel_options)
    revalidation,oracle=_final_bind(solution,terminal,options,production_options,source_mode)
    context=_final_context(solution,terminal,production_options)
    issues=FinalObservationIssue[]
    views=_final_views(solution,context,issues)
    positivity=try
        positivity_diagnostics(solution;integrate=true)
    catch error
        _final_expected_error(error) || rethrow()
        _final_issue!(issues,:positivity,:positivity_measurement_failed,sprint(showerror,error))
    end
    assessment=solution_scientific_assessment(solution)
    # Only this existing getter field aliases its source. Fresh getter arrays transfer ownership.
    assessment["final_fixed_point_audit"]=deepcopy(assessment["final_fixed_point_audit"])
    summary=FinalScientificSummary(solution.status,solution_quality(solution),solution.converged,
        assessment["registry_version"],assessment["iterative_converged"],assessment["fixed_hartree_converged"],
        assessment["physical_gates_passed"],assessment["stationary_candidate_accepted"],
        assessment["scientific_accepted"],Symbol(assessment["discretization_verified"]),
        assessment["scientific_acceptance_reason"],assessment["outer_confirmation_count"])
    return FinalObservationSnapshot(source_mode,terminal.identity,terminal.integrity,context,views,
        solution_metric_diagnostics(solution),positivity,summary,assessment,revalidation,oracle,Tuple(issues))
end
