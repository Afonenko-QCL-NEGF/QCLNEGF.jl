# Serialization-neutral, scoped in-memory snapshots. Contracts owns file envelopes.
const ExactDigest = NTuple{32,UInt8}
const EXACT_STATE_VERSION = UInt16(1)

"""Explicit numerical runtime facts; unknown facts cannot qualify an exact restart."""
struct ExactExecutionIdentity
    core_code::String
    julia::String
    blas::String
    fftw::String
    scalar_type::String
    index_type::String
    reduction::String
    worker_count::Int
    blas_threads::Int
    energy_chunk::Int
    hilbert_columns::Int
    residual_chunk::Int
    parallel_backend::Symbol
    function ExactExecutionIdentity(; core_code,julia,blas,fftw,scalar_type,index_type,reduction,
        worker_count,blas_threads,energy_chunk,hilbert_columns,residual_chunk,parallel_backend)
        strings=String.((core_code,julia,blas,fftw,scalar_type,index_type,reduction))
        all(x->!isempty(strip(x)) && lowercase(strip(x)) ∉ ("unknown","not_measured"),strings) ||
            throw(ArgumentError("execution identity has missing source/runtime/type/reduction facts"))
        counts=Int.((worker_count,blas_threads,energy_chunk,hilbert_columns,residual_chunk))
        all(>(0),counts) || throw(ArgumentError("execution widths/chunks must be positive"))
        parallel_backend in (:threads,:blas) || throw(ArgumentError("execution backend is unknown"))
        new(strings...,counts...,parallel_backend)
    end
end

struct ExactRestartIdentity
    contract_version::UInt16
    problem::ExactDigest
    policy::ExactDigest
    execution::ExactDigest
end
Base.:(==)(a::ExactRestartIdentity,b::ExactRestartIdentity) =
    all(name->getfield(a,name)==getfield(b,name),fieldnames(ExactRestartIdentity))

"""Finite value grammar for Core warnings/audit records, not a persistence protocol."""
abstract type CoreDatum end
struct CoreScalar <: CoreDatum
    value::Union{Nothing,Bool,Int64,Float64,ComplexF64,String,Symbol}
end
struct CoreArray{T<:Union{Bool,Int64,Float64,ComplexF64}} <: CoreDatum
    shape::Tuple{Vararg{Int}}
    values::Vector{T}
    function CoreArray(array::AbstractArray{T}) where {T<:Union{Bool,Int64,Float64,ComplexF64}}
        new{T}(size(array),copy(vec(array)))
    end
end
struct CoreList <: CoreDatum
    values::Tuple{Vararg{CoreDatum}}
end
struct CoreRecord <: CoreDatum
    entries::Tuple{Vararg{Pair{String,CoreDatum}}}
    function CoreRecord(entries::Tuple)
        keys=String[String(first(entry)) for entry in entries]
        length(unique(keys))==length(keys) || throw(ArgumentError("record has duplicate keys"))
        pairs=Pair{String,CoreDatum}[
            String(first(entry))=>_core_datum(last(entry)) for entry in entries]
        sort!(pairs;by=first)
        new(Tuple(pairs))
    end
end
CoreRecord(value::AbstractDict)=CoreRecord(Tuple(pairs(value)))
CoreRecord(value::NamedTuple)=CoreRecord(Tuple(pairs(value)))
function Base.getindex(record::CoreRecord,key::Union{String,Symbol})
    for pair in record.entries
        first(pair)==String(key) && return last(pair)
    end
    throw(KeyError(key))
end
Base.haskey(record::CoreRecord,key::Union{String,Symbol}) =
    any(pair->first(pair)==String(key),record.entries)
_core_datum(value::CoreDatum)=deepcopy(value)
for atom in (Nothing,Bool,Int64,Float64,ComplexF64,String,Symbol)
    @eval _core_datum(value::$atom)=CoreScalar(value)
end
_core_datum(value::Integer)=CoreScalar(Int64(value))
_core_datum(value::AbstractDict)=CoreRecord(value)
_core_datum(value::NamedTuple)=CoreRecord(value)
_core_datum(value::AbstractArray{T}) where {T<:Union{Bool,Int64,Float64,ComplexF64}}=CoreArray(value)
_core_datum(value::AbstractVector{T}) where {T<:Union{Bool,Int64,Float64,ComplexF64}}=CoreArray(value)
_core_datum(value::Union{AbstractVector,Tuple})=CoreList(Tuple(_core_datum(x) for x in value))
_core_datum(value)=throw(ArgumentError("unsupported Core evidence value $(typeof(value))"))
_core_plain(value::CoreScalar)=value.value
_core_plain(value::CoreArray)=reshape(copy(value.values),value.shape)
_core_plain(value::CoreList)=[_core_plain(x) for x in value.values]
_core_plain(value::CoreRecord)=Dict{String,Any}(first(x)=>_core_plain(last(x)) for x in value.entries)

struct SolverProgress
    completed_outer::Int
    current_outer::Int
    accepted_scba::Int
    outer_strict_streak::Int
    inner_strict_streak::Int
    diagnostic_streak::Int
end
SolverProgress(;completed_outer,current_outer,accepted_scba,outer_strict_streak,
    inner_strict_streak,diagnostic_streak)=SolverProgress(completed_outer,current_outer,
    accepted_scba,outer_strict_streak,inner_strict_streak,diagnostic_streak)

struct SCBAContinuation
    decision::SCBAStopDecision
    candidate::Union{Nothing,FixedPointAuditCandidate}
    next_previous_current::Float64
    next_previous_populations::Vector{Float64}
    seed_mu_scaled::Union{Nothing,Float64}
    seed_number_ratio::Union{Nothing,Float64}
    witness_retention::NTuple{3,Int}
end
function SCBAContinuation(;decision,candidate,next_previous_current,next_previous_populations,
    seed_mu_scaled,seed_number_ratio,witness_retention)
    !decision.terminal && candidate===nothing &&
        throw(ArgumentError("continuation.candidate is missing for pending SCBA mixing"))
    SCBAContinuation(decision,deepcopy(candidate),Float64(next_previous_current),
        copy(next_previous_populations),seed_mu_scaled,seed_number_ratio,witness_retention)
end
struct OuterContinuation
    next_U::Vector{Float64}
    next_action::Symbol
end
function OuterContinuation(;next_U,next_action)
    next_action in (:next_outer,:finalize) || throw(ArgumentError("continuation.next_action is unknown"))
    OuterContinuation(copy(next_U),next_action)
end

struct AcceptanceContext
    backend::Symbol
    algorithms::AlgorithmOptions
    hilbert_roundoff_measured::Union{Nothing,Bool}
    map_kernels_shared_with_solver::Union{Nothing,Bool}
    final_state_alignment::Symbol
end
AcceptanceContext(;backend,algorithms,hilbert_roundoff_measured,map_kernels_shared_with_solver,
    final_state_alignment)=AcceptanceContext(backend,algorithms,hilbert_roundoff_measured,
    map_kernels_shared_with_solver,final_state_alignment)

struct ProblemRebuildRecipe
    physical::PhysicalParameters
    numerical::NumericalParameters
    scales::ScaleSystem
    scattering::ScatteringOptions
    models::PhysicalModelOptions
    algorithms::AlgorithmOptions
    kernel_options::ProductionKernelOptions
end
ProblemRebuildRecipe(;physical,numerical,scales,scattering,models,algorithms,kernel_options)=
    ProblemRebuildRecipe(deepcopy(physical),numerical,scales,scattering,models,algorithms,kernel_options)
"""Required evidence of a domain skipped by the existing representation preflight."""
struct DomainPreflightEvidence
    discretization_id::String
    basis_id::String
    E_min_eV::Float64
    E_max_eV::Float64
    NE::Int
    energy_min_scaled::Float64
    energy_max_scaled::Float64
    required_min_scaled::Float64
    required_max_scaled::Float64
    reason::String
    seed_transfer::String
end
DomainPreflightEvidence(;discretization_id,basis_id,E_min_eV,E_max_eV,NE,
    energy_min_scaled,energy_max_scaled,required_min_scaled,required_max_scaled,reason,seed_transfer)=
    DomainPreflightEvidence(discretization_id,basis_id,E_min_eV,E_max_eV,NE,
        energy_min_scaled,energy_max_scaled,required_min_scaled,required_max_scaled,reason,seed_transfer)

"""Required evidence of a measured completed domain; IEEE failed diagnostics are retained."""
struct DomainMeasuredEvidence
    discretization_id::String
    basis_id::String
    hartree_id::String
    E_min_eV::Float64
    E_max_eV::Float64
    NE::Int
    energy_min_scaled::Float64
    energy_max_scaled::Float64
    required_min_scaled::Float64
    required_max_scaled::Float64
    tail_low::Float64
    tail_high::Float64
    edge_gamma::Float64
    edge_spectral::Float64
    scba_quality::Symbol
    seconds::Float64
    status::String
    strict_accepted::Bool
    metrics::Dict{Symbol,Float64}
    current_A_per_m2::Float64
    momentum_tail_unresolved::Bool
    seed_transfer::String
    reason::String
end
DomainMeasuredEvidence(;discretization_id,basis_id,hartree_id,E_min_eV,E_max_eV,NE,
    energy_min_scaled,energy_max_scaled,required_min_scaled,required_max_scaled,
    tail_low,tail_high,edge_gamma,edge_spectral,scba_quality,seconds,
    status,strict_accepted,metrics,current_A_per_m2,momentum_tail_unresolved,seed_transfer,reason)=
    DomainMeasuredEvidence(discretization_id,basis_id,hartree_id,E_min_eV,E_max_eV,NE,
        energy_min_scaled,energy_max_scaled,required_min_scaled,required_max_scaled,
        tail_low,tail_high,edge_gamma,edge_spectral,scba_quality,seconds,
        status,strict_accepted,Dict{Symbol,Float64}(Symbol(k)=>v for (k,v) in metrics),
        current_A_per_m2,momentum_tail_unresolved,seed_transfer,reason)

struct DomainRevision
    revision::Int
    source::Symbol
    evidence::Union{DomainPreflightEvidence,DomainMeasuredEvidence}
end
function DomainRevision(;revision,source,evidence)
    ((source===:preflight_estimate && evidence isa DomainPreflightEvidence) ||
     (source===:measured_state && evidence isa DomainMeasuredEvidence)) ||
        throw(ArgumentError("adaptation.revision.source differs from its typed evidence"))
    DomainRevision(revision,source,deepcopy(evidence))
end
struct DomainAdaptationState
    policy::DomainAdaptationPolicy
    kernel_options::ProductionKernelOptions
    expansions::Int
    elapsed_domain_seconds::Float64
    current_problem::ExactDigest
    current_problem_recipe::ProblemRebuildRecipe
    revisions::Tuple{Vararg{DomainRevision}}
    next_action::Symbol
    decision_id::Int
    decision_committed::Bool
    next_problem_recipe::Union{Nothing,ProblemRebuildRecipe}
end
function DomainAdaptationState(;policy,kernel_options,expansions,elapsed_domain_seconds,current_problem,current_problem_recipe,revisions,
    next_action,decision_id,decision_committed,next_problem_recipe)
    next_action in (:assess_measured_domain,:cold_rebuild,:finished) ||
        throw(ArgumentError("adaptation.next_action is unknown"))
    next_action===:cold_rebuild && next_problem_recipe===nothing &&
        throw(ArgumentError("adaptation.next_problem_recipe is missing"))
    DomainAdaptationState(policy,kernel_options,expansions,elapsed_domain_seconds,current_problem,deepcopy(current_problem_recipe),deepcopy(Tuple(revisions)),
        next_action,decision_id,decision_committed,deepcopy(next_problem_recipe))
end

struct FinalEvidenceSnapshot
    report::ConvergenceReport
    original_converged::Bool
    original_status::Symbol
    observations::CoreRecord
end
FinalEvidenceSnapshot(;report,original_converged,original_status,observations)=
    FinalEvidenceSnapshot(deepcopy(report),original_converged,original_status,deepcopy(observations))

struct ScientificHistorySnapshot
    action::Symbol
    stage::Symbol
    domain_revision::Int
    discretization_id::String
    basis_id::String
    outer_iteration::Int
    inner_iteration::Union{Nothing,Int}
    sequence::Int
    converged::Bool
    status::Symbol
    quality::Symbol
    strict_streak::Int
    diagnostic_streak::Int
    units::Tuple{Vararg{Pair{Symbol,String}}}
    missing_reasons::Tuple{Vararg{Pair{Symbol,Symbol}}}
    row::Union{SCBAIteration,OuterIteration}
end
function _scientific_row_metadata(row,status)
    units=Pair{Symbol,String}[]
    missing=Pair{Symbol,Symbol}[]
    first_iteration=row isa SCBAIteration ? row.ν==1 : row.μ==1
    for name in fieldnames(typeof(row))
        value=getfield(row,name)
        if value isa Real || value isa Vector{Float64}
            unit=name===:J ? "A/m^2" : name===:density ? "L0^3 n" :
                name===:populations ? "target_sheet_fraction" : name===:r_U ? "E0" :
                name in (:raw_charge,:target_charge,:normalized_charge) ? "L0^2 n_sheet" : "1"
            push!(units,name=>unit)
            finite=value isa Real ? isfinite(value) : all(isfinite,value)
            if !finite
                previous=first_iteration && name in (:r_Jchange,:r_population,:lambda_change,:r_n)
                push!(missing,name=>(previous ? :no_previous_iteration : :non_finite))
            end
        elseif value===nothing && name in (:witness,:physical_markers)
            failed=status in (:invalid_candidate,:nonfinite,:nonfinite_metrics)
            push!(missing,name=>(failed ? :measurement_failed :
                name===:physical_markers ? :not_sampled : :not_measured))
        end
    end
    if row isa SCBAIteration && row.witness!==nothing &&
       isempty(row.witness.matrices) && isempty(row.witness.eigenvector)
        push!(missing,:witness_payload=>:not_retained_by_policy)
    end
    return Tuple(units),Tuple(missing)
end
function ScientificHistorySnapshot(;stage,domain_revision,discretization_id,basis_id,
    outer_iteration,inner_iteration,sequence,
    converged,status,quality,strict_streak,diagnostic_streak,row,action::Symbol=:append)
    action in (:append,:witness_update) || throw(ArgumentError("history.action is unknown"))
    ((stage===:scba && row isa SCBAIteration) || (stage===:outer && row isa OuterIteration)) ||
        throw(ArgumentError("history.stage does not match row type"))
    domain_revision>=0 && !isempty(discretization_id) && !isempty(basis_id) ||
        throw(ArgumentError("history domain context is missing or invalid"))
    units,missing=_scientific_row_metadata(row,status)
    ScientificHistorySnapshot(action,stage,domain_revision,discretization_id,basis_id,
        outer_iteration,inner_iteration,sequence,converged,status,
        quality,strict_streak,diagnostic_streak,units,missing,deepcopy(row))
end

struct ExactSolverState
    identity::ExactRestartIdentity
    boundary::Symbol
    published_U::Vector{Float64}
    published_n::Vector{Float64}
    scba::SCBAResult
    mixer_position::MixerPosition
    family_order::Tuple{Vararg{Symbol}}
    outer_history::Vector{OuterIteration}
    progress::SolverProgress
    continuation::Union{Nothing,SCBAContinuation,OuterContinuation}
    status::Symbol
    warnings::Tuple{Vararg{CoreRecord}}
    policy_evaluations::Tuple{Vararg{CoreRecord}}
    adaptation::Union{Nothing,DomainAdaptationState}
    measurement_context::AcceptanceContext
    final_evidence::Union{Nothing,FinalEvidenceSnapshot}
    integrity::ExactDigest
end
function ExactSolverState(;identity,boundary,published_U,published_n,scba,mixer_position,family_order,outer_history,
    progress,continuation,status,warnings,policy_evaluations,adaptation,measurement_context,final_evidence,
    integrity=nothing)
    boundary in (:scba_accepted,:outer_accepted,:domain_final,:final) ||
        throw(ArgumentError("state.boundary is unknown"))
    payload=deepcopy((identity,boundary,published_U,published_n,scba,mixer_position,Tuple(family_order),outer_history,progress,
        continuation,status,Tuple(warnings),Tuple(policy_evaluations),adaptation,measurement_context,final_evidence))
    calculated=_exact_digest(payload)
    integrity===nothing || integrity==calculated || throw(ArgumentError("state.integrity differs"))
    ExactSolverState(payload...,calculated)
end

"""Opaque owning handle. Construct with restore_exact_state; do not inspect its fields."""
struct RestoredSolverState
    problem::NEGFProblem
    state::ExactSolverState
    options::SolverOptions
    production_options::ProductionOptions
    execution::ExactExecutionIdentity
    domain_adaptation::DomainAdaptationPolicy
    kernel_options::Union{Nothing,ProductionKernelOptions}
    function RestoredSolverState(problem,state,options,production_options,execution,
                                domain_adaptation,kernel_options,::Val{:validated})
        new(problem,deepcopy(state),options,production_options,execution,domain_adaptation,kernel_options)
    end
end
