"""Checked numeric-payload/metadata allowances, excluding compiler/allocator RSS."""
struct FinalObservationMemoryEstimate
    dimensions::NTuple{4,Int}
    caller_resident_bytes::Int
    context_bytes::Int
    view_bytes::Int
    positivity_bytes::Int
    metadata_bytes::Int
    workspace_bytes::Int
    existing_validation_bytes::Int
    peak_additional_bytes::Int
end
struct FinalObservationScales
    E0_eV::Float64
    L0_m::Float64
    J0_A_per_m2::Float64
end
struct FinalObservationInputs
    F_bias_V_per_m::Float64
    V_period_V::Float64
    T_L_K::Float64
    T_LO_K::Float64
    E_ref_eV::Float64
    N_dop_2D_per_m2::Float64
    f_ion::Float64
    z0_m::Float64
    period_length_m::Float64
    g_s::Int
end
"""Detached actual-domain context. Construction flags are not original Runner requests."""
struct FinalObservationContext
    dimensions::NTuple{4,Int}
    domain_revision::Int
    discretization_id::String
    basis_id::String
    x_scaled::Vector{Float64}
    wX::Vector{Float64}
    energy_scaled::Vector{Float64}
    wE::Vector{Float64}
    momentum_scaled::Vector{Float64}
    wk::Vector{Float64}
    z_nm::Vector{Float64}
    energy_eV::Vector{Float64}
    k_per_nm::Vector{Float64}
    state_index::Vector{Int64}
    trusted_energy::BitVector
    scales::FinalObservationScales
    inputs::FinalObservationInputs
    problem_scattering::ScatteringOptions
    enabled_scattering::Tuple{Vararg{Symbol}}
    selfenergy_families::Tuple{Vararg{Symbol}}
    physical_models::PhysicalModelOptions
    algorithms::AlgorithmOptions
    energy_shift_discretization::Symbol
    localization::Symbol
    boundary_condition::Symbol
    basis_contract::Symbol
    energy_reference::Symbol
    effective_gauge::Symbol
    function FinalObservationContext(dimensions::NTuple{4,Int},
        domain_revision::Int,
        discretization_id::String,
        basis_id::String,
        x_scaled::Vector{Float64},
        wX::Vector{Float64},
        energy_scaled::Vector{Float64},
        wE::Vector{Float64},
        momentum_scaled::Vector{Float64},
        wk::Vector{Float64},
        z_nm::Vector{Float64},
        energy_eV::Vector{Float64},
        k_per_nm::Vector{Float64},
        state_index::Vector{Int64},
        trusted_energy::BitVector,
        scales::FinalObservationScales,
        inputs::FinalObservationInputs,
        problem_scattering::ScatteringOptions,
        enabled_scattering::Tuple{Vararg{Symbol}},
        selfenergy_families::Tuple{Vararg{Symbol}},
        physical_models::PhysicalModelOptions,
        algorithms::AlgorithmOptions,
        energy_shift_discretization::Symbol,
        localization::Symbol,
        boundary_condition::Symbol,
        basis_contract::Symbol,
        energy_reference::Symbol,
        effective_gauge::Symbol)
        _final_check_context(;dimensions,domain_revision,discretization_id,basis_id,x_scaled,wX,energy_scaled,
            wE,momentum_scaled,wk,z_nm,energy_eV,k_per_nm,state_index,trusted_energy,scales,inputs,
            problem_scattering,enabled_scattering,selfenergy_families,physical_models,algorithms,
            energy_shift_discretization,localization,boundary_condition,basis_contract,energy_reference,
            effective_gauge)
        new(dimensions,domain_revision,discretization_id,basis_id,x_scaled,wX,energy_scaled,wE,momentum_scaled,
            wk,z_nm,energy_eV,k_per_nm,state_index,trusted_energy,scales,inputs,problem_scattering,
            enabled_scattering,selfenergy_families,physical_models,algorithms,energy_shift_discretization,
            localization,boundary_condition,basis_contract,energy_reference,effective_gauge)
    end
end
"""Owned stored views with named physical units; `nothing` requires a matching snapshot issue."""
struct StoredFinalObservations
    potential_structure_eV::Vector{Float64}
    potential_external_eV::Vector{Float64}
    potential_hartree_eV::Vector{Float64}
    potential_total_eV::Vector{Float64}
    density_per_m3::Vector{Float64}
    spatial_energy_density_per_eV_m3::Union{Nothing,Matrix{Float64}}
    spectral_trace_per_eV::Matrix{Float64}
    occupied_trace_per_eV::Matrix{Float64}
    localized_levels_eV::Union{Nothing,Vector{Float64}}
    effective_levels_eV::Union{Nothing,Vector{Float64}}
    localized_wavefunctions_per_sqrt_nm::Matrix{ComplexF64}
    effective_wavefunctions_per_sqrt_nm::Union{Nothing,Matrix{ComplexF64}}
    effective_transform::Union{Nothing,Matrix{ComplexF64}}
    sheet_density_matrix_per_m2::Matrix{ComplexF64}
    effective_density_matrix_per_m2::Union{Nothing,Matrix{ComplexF64}}
    current_density_A_per_m2::Float64
    sheet_density_per_m2::Float64
    function StoredFinalObservations(potential_structure_eV::Vector{Float64},
        potential_external_eV::Vector{Float64},
        potential_hartree_eV::Vector{Float64},
        potential_total_eV::Vector{Float64},
        density_per_m3::Vector{Float64},
        spatial_energy_density_per_eV_m3::Union{Nothing,Matrix{Float64}},
        spectral_trace_per_eV::Matrix{Float64},
        occupied_trace_per_eV::Matrix{Float64},
        localized_levels_eV::Union{Nothing,Vector{Float64}},
        effective_levels_eV::Union{Nothing,Vector{Float64}},
        localized_wavefunctions_per_sqrt_nm::Matrix{ComplexF64},
        effective_wavefunctions_per_sqrt_nm::Union{Nothing,Matrix{ComplexF64}},
        effective_transform::Union{Nothing,Matrix{ComplexF64}},
        sheet_density_matrix_per_m2::Matrix{ComplexF64},
        effective_density_matrix_per_m2::Union{Nothing,Matrix{ComplexF64}},
        current_density_A_per_m2::Float64,
        sheet_density_per_m2::Float64)
        _final_check_views(;potential_structure_eV,potential_external_eV,potential_hartree_eV,
            potential_total_eV,density_per_m3,spatial_energy_density_per_eV_m3,spectral_trace_per_eV,
            occupied_trace_per_eV,localized_levels_eV,effective_levels_eV,localized_wavefunctions_per_sqrt_nm,
            effective_wavefunctions_per_sqrt_nm,effective_transform,sheet_density_matrix_per_m2,
            effective_density_matrix_per_m2,current_density_A_per_m2,sheet_density_per_m2)
        new(potential_structure_eV,potential_external_eV,potential_hartree_eV,potential_total_eV,density_per_m3,
            spatial_energy_density_per_eV_m3,spectral_trace_per_eV,occupied_trace_per_eV,localized_levels_eV,
            effective_levels_eV,localized_wavefunctions_per_sqrt_nm,effective_wavefunctions_per_sqrt_nm,
            effective_transform,sheet_density_matrix_per_m2,effective_density_matrix_per_m2,
            current_density_A_per_m2,sheet_density_per_m2)
    end
end
struct FinalScientificSummary
    status::Symbol
    quality::Symbol
    converged::Bool
    registry_version::String
    iterative_converged::Bool
    fixed_hartree_converged::Bool
    physical_gates_passed::Bool
    stationary_candidate_accepted::Bool
    scientific_accepted::Union{Nothing,Bool}
    discretization_verified::Symbol
    scientific_acceptance_reason::String
    outer_confirmation_count::Int
end
struct FinalObservationIssue
    target::Symbol
    status::Symbol
    reason::Symbol
    detail::String
    function FinalObservationIssue(target::Symbol,status::Symbol,reason::Symbol,detail::AbstractString)
        status in (:error,:not_measured) || throw(ArgumentError("observation issue status is unknown"))
        new(target,status,reason,String(detail))
    end
end
"""Owning result of final_observation_snapshot, outside the exact-state persistence graph.

Array fields may be mutated/retained by the recipient. Positional typed values
alone attest neither a producer nor scientific validity; the capture port checks
source shapes, terminal identity and raw-bit binding before producing this value.
"""
struct FinalObservationSnapshot
    source_mode::Symbol
    identity::ExactRestartIdentity
    terminal_integrity::ExactDigest
    context::FinalObservationContext
    views::StoredFinalObservations
    metrics::SolutionMetricDiagnostics
    positivity::Union{Nothing,PositivityDiagnostics}
    summary::FinalScientificSummary
    assessment::Dict{String,Any}
    revalidation::Union{Nothing,RevalidationReceipt}
    independent_oracle::Union{Nothing,IndependentOracleReceipt}
    issues::Tuple{Vararg{FinalObservationIssue}}
    function FinalObservationSnapshot(source_mode::Symbol,
        identity::ExactRestartIdentity,
        terminal_integrity::ExactDigest,
        context::FinalObservationContext,
        views::StoredFinalObservations,
        metrics::SolutionMetricDiagnostics,
        positivity::Union{Nothing,PositivityDiagnostics},
        summary::FinalScientificSummary,
        assessment::Dict{String,Any},
        revalidation::Union{Nothing,RevalidationReceipt},
        independent_oracle::Union{Nothing,IndependentOracleReceipt},
        issues::Tuple{Vararg{FinalObservationIssue}})
        _final_check_snapshot(;source_mode,identity,terminal_integrity,context,views,metrics,positivity,summary,
            assessment,revalidation,independent_oracle,issues)
        new(source_mode,identity,terminal_integrity,context,views,metrics,positivity,summary,assessment,
            revalidation,independent_oracle,issues)
    end
end
