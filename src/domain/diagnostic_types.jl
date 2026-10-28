"""Detached summary values and declared strict/optional reporting-band limits.

These metadata do not assert measurement, convergence or scientific acceptance.
"""
struct SolutionMetricDiagnostics
    summary_metrics::Dict{Symbol,Float64}
    strict_limits::Dict{Symbol,Float64}
    approximate_limits::Union{Nothing,Dict{Symbol,Float64}}
end

"""Immutable worst-block PSD evidence in scaled solver units.

`status=:error` preserves failed IEEE diagnostics, including legacy compatibility
defaults; none of its scalar fields may be interpreted as a measured pass.
"""
struct PositivityWitness
    status::Symbol
    matrix_kind::Symbol
    energy_index::Int
    momentum_index::Int
    minimum_eigenvalue::Float64
    block_norm::Float64
    floor::Float64
    ratio::Float64
    absolute_defect::Float64
    relative_defect::Float64
    backward_error::Float64
    hermiticity_defect::Float64
    construction_scale::Float64
end

"""Complete-grid matrix integrals, negative weights and invalid-block counts.

Spectral integrals/eigenvalues are dimensionless. The negative energy integral
has `(kind,k)` axes in scaled solver units (broadening carries energy squared).
The five phase-space weights have m⁻² units and omit the broadening kind.
All arrays belong to this result, independently of the observed solution.
"""
struct MatrixIntegralDiagnostics
    spectral_integral::Array{ComplexF64,3}
    spectral_integral_eigenvalues::Matrix{Float64}
    spectral_trace_cumulative::Matrix{Float64}
    negative_energy_integral_scaled::Matrix{Float64}
    negative_block_count::NTuple{6,Int}
    beyond_backward_error_block_count::NTuple{6,Int}
    invalid_block_count::NTuple{6,Int}
    total_blocks_per_kind::Int
    negative_phase_space_weight_per_m2::NTuple{5,Float64}
end

"""Six-kind PSD observations; local threshold flags are not global acceptance.

`integrals=nothing` means not requested. `scan_order` declares the strict-`>`
tie order of the selected witnesses, rather than silently changing coordinates.
"""
struct PositivityDiagnostics
    metric_version::Int
    matrix_kinds::NTuple{6,Symbol}
    witnesses::NTuple{6,PositivityWitness}
    strict_ratio_limit::Float64
    strict_absolute_limits_scaled::NTuple{6,Float64}
    strict_threshold_satisfied::NTuple{6,Bool}
    scan_order::Symbol
    integrals::Union{Nothing,MatrixIntegralDiagnostics}
end

"""Measured fixed-basis subset of static validation, with detached limits."""
struct BasisValidationDiagnostics
    report::ConvergenceReport
    limits::Dict{Symbol,Float64}
end

"""Detached projected eigenpairs and the unchanged effective-level gate.

Eigenvalues use the input Hamiltonian's energy coordinates; they are diagnostic
one-particle levels, not interacting spectral peaks.
"""
struct ProjectedHamiltonianSpectrum
    eigenvalues::Vector{Float64}
    eigenvectors::Matrix{ComplexF64}
    hermiticity_residual::Float64
    eigen_residual::Float64
    norm_floor::Float64
    hermiticity_limit::Float64
end

"""Owned optimized resolvent `(E,k,a,b)` and condition/scale arrays `(E,k)`."""
struct ProductionRetardedGreenEvaluation
    G_R::Array{ComplexF64,4}
    condition_number::Matrix{Float64}
    dyson_scale::Matrix{Float64}
end

"""Owned LO lesser/greater collision self-energies `(E,k,a,b)`, in E₀ units."""
struct LOCollisionEvaluation
    lesser::Array{ComplexF64,4}
    greater::Array{ComplexF64,4}
end

"""Sampled modal Γ/ΔE observations with unchanged status/reason and NaN semantics."""
struct SampledLinewidthDiagnostics
    q10::Float64
    q50::Float64
    q90::Float64
    underresolved::Float64
    blocks::Int
    status::Symbol
    sampling_method::Symbol
end
