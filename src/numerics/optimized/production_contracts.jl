"""
    EnergyShiftPlan

Sparse, non-periodic representation of the piecewise-linear map
`X(ε) -> X(ε + δ)`.  Every valid target energy uses at most two source
nodes.  Invalid rows are exactly zero; energy wrap-around is impossible.

The plan is algebraically equivalent to [`build_shift_matrix`](@ref) but
requires `O(N_E)` rather than `O(N_E^2)` storage and work.

See [Field-periodic closure](@ref theory-periodicity) and
[the production sparse-shift construction](@ref theory-production).
"""
struct EnergyShiftPlan <: AbstractMatrix{Float64}
    left::Vector{Int}
    right::Vector{Int}
    weight_right::Vector{Float64}
    valid::BitVector
    δ::Float64
end

"""One serialization-neutral scalar attached to a [`SolverEvent`](@ref)."""
struct SolverMetric
    name::Symbol
    value::Union{Float64,Int64,Bool,String}
    unit::String
end

# The generated union-field constructor and the abstract convenience methods
# below are otherwise ambiguous for, for example, `(Float64, String)`: one is
# more specific in the value argument and the other in the unit argument.
# Exact forwarding methods make every supported scalar type unambiguous while
# preserving the serialization-neutral event contract.
for T in (Float64, Int64, Bool, String)
    @eval SolverMetric(name::Symbol, value::$T, unit::String) = invoke(
        SolverMetric,
        Tuple{Symbol,Union{Float64,Int64,Bool,String},String},
        name,
        value,
        unit,
    )
end

SolverMetric(name::Symbol, value::AbstractFloat, unit::AbstractString = "") =
    SolverMetric(name, Float64(value), String(unit))
SolverMetric(name::Symbol, value::Integer, unit::AbstractString = "") =
    SolverMetric(name, Int64(value), String(unit))
SolverMetric(name::Symbol, value::Bool, unit::AbstractString = "") =
    SolverMetric(name, value, String(unit))
SolverMetric(name::Symbol, value::AbstractString, unit::AbstractString = "") =
    SolverMetric(name, String(value), String(unit))

"""Neutral numerical event consumed by an application-supplied callback."""
struct SolverEvent
    action::Symbol
    stage::Symbol
    label::String
    status::Symbol
    iteration::Union{Nothing,Int}
    total::Union{Nothing,Int}
    metrics::Vector{SolverMetric}
    message::String
end

_solver_metrics(::Nothing) = SolverMetric[]
_solver_metrics(metric::SolverMetric) = SolverMetric[metric]
_solver_metrics(metrics) = SolverMetric[metric for metric in metrics]

"""
    ProductionOptions

Controls memory-bounded production operators, checkpoint cadence, and neutral
progress events without changing the physical model. `parallel_backend=:threads`
partitions independent energy chunks and matrix blocks among Julia workers;
use it with multi-threaded Julia and one BLAS thread. `:blas` uses one outer
worker for BLAS/LAPACK-heavy contraction, Dyson/Keldysh, embedding, and
residual blocks, while BLAS may use its own threads; private-buffer FFT and
lightweight elementwise phases may still use Julia threads. `worker_count=0`
lets the application resource planner use all available Julia threads, while a
positive value caps the logical worker count. `residual_chunk` sets
deterministic residual task granularity, and `phase_timing=true` emits
per-iteration timing telemetry.

See [Production backend](@ref theory-production),
[numerical cost and memory](@ref theory-cost), and the
[state persistence contract](https://github.com/AfonenkoA/QCLNEGFRunner.jl/blob/main/docs/src/user/results.md).
"""
Base.@kwdef struct ProductionOptions
    # The application resource planner replaces this sentinel before every
    # production run.  Low-level numerics never probes the host or cgroup.
    memory_budget_bytes::Int = typemax(Int)
    energy_chunk::Int = 64
    hilbert_columns::Int = 32
    parallel_backend::Symbol = :threads
    worker_count::Int = 0
    residual_chunk::Int = 256
    phase_timing::Bool = false
    physics_markers::SCBAPhysicsMarkerPolicy = SCBAPhysicsMarkerPolicy()
    verify_fft_roundoff::Bool = true
    checkpoint_every_outer::Int = 5
    checkpoint_every_scba::Int = 25
    progress_every_scba::Int = 5
    progress_every_outer::Int = 1
    algorithms::AlgorithmOptions = AlgorithmOptions()
    event_sink::Union{Nothing,Function} = nothing
    checkpoint_request::Union{Nothing,Function} = nothing
    phase_request::Union{Nothing,Function} = nothing
    outer_iteration::Int = 0
end

"""Return a copy of `options` with explicitly named fields replaced.

This is the single copy/update constructor used by configuration, studies,
auto-tuning, and progress wiring. Keeping it here prevents positional copies
of the production policy from drifting when a field is added.
"""
function with_production_options(
    options::ProductionOptions;
    memory_budget_bytes = options.memory_budget_bytes,
    energy_chunk = options.energy_chunk,
    hilbert_columns = options.hilbert_columns,
    parallel_backend = options.parallel_backend,
    worker_count = options.worker_count,
    residual_chunk = options.residual_chunk,
    phase_timing = options.phase_timing,
    physics_markers = options.physics_markers,
    verify_fft_roundoff = options.verify_fft_roundoff,
    checkpoint_every_outer = options.checkpoint_every_outer,
    checkpoint_every_scba = options.checkpoint_every_scba,
    progress_every_scba = options.progress_every_scba,
    progress_every_outer = options.progress_every_outer,
    algorithms = options.algorithms,
    event_sink = options.event_sink,
    checkpoint_request = options.checkpoint_request,
    phase_request = options.phase_request,
    outer_iteration = options.outer_iteration,
)
    return ProductionOptions(;
        memory_budget_bytes,
        energy_chunk,
        hilbert_columns,
        parallel_backend,
        worker_count,
        residual_chunk,
        phase_timing,
        physics_markers,
        verify_fft_roundoff,
        checkpoint_every_outer,
        checkpoint_every_scba,
        progress_every_scba,
        progress_every_outer,
        algorithms,
        event_sink,
        checkpoint_request,
        phase_request,
        outer_iteration,
    )
end

"""
Conservative resident/peak-memory and work estimate for one solve.

See [Numerical scales and cost](@ref theory-cost) and the
[production RAM contract](@ref theory-production).
"""
struct ProductionMemoryEstimate
    resident_bytes::Int
    peak_bytes::Int
    memory_budget_bytes::Int
    kernel_flops_per_candidate::Float64
    breakdown::Dict{Symbol,Int}
end

function Base.show(io::IO, estimate::ProductionMemoryEstimate)
    gib = 1024.0^3
    print(
        io,
        "ProductionMemoryEstimate(resident=",
        round(estimate.resident_bytes / gib; digits = 2),
        " GiB, peak=",
        round(estimate.peak_bytes / gib; digits = 2),
        " GiB, limit=",
        round(estimate.memory_budget_bytes / gib; digits = 2),
        " GiB, ",
        "kernel_work=",
        round(estimate.kernel_flops_per_candidate; sigdigits = 3),
        " real-FLOP-equivalents)",
    )
end

"""
Common dispatch interface for production scattering operators.

Concrete implementations represent the same flattened six-index contraction
literally, as a dense matrix, through an exact momentum-independent factor,
or through a controlled truncated SVD.  The selected representation and its
scientific impact are recorded by [`algorithm_manifest`](@ref QCLNEGF.algorithm_manifest).

See [Production backend](@ref theory-production),
[Microscopic kernels](@ref theory-kernels), and
[Optimization decision tree](@ref optimization-decision-tree).
"""
abstract type AbstractProductionKernel end

"""
Marker for the reference six-index contraction.  The tensor is not duplicated:
[`ProductionCache`](@ref) retains the normalized source array and the literal
loops consume it directly.
"""
struct LiteralProductionKernel <: AbstractProductionKernel
    N_k::Int
    N_b::Int
end

"""
Dense exact scattering map from `(k′,c,d)` to `(k,a,b)`.

See [Microscopic kernels](@ref theory-kernels) and the
[production BLAS contraction](@ref theory-production).
"""
struct DenseProductionKernel <: AbstractProductionKernel
    matrix::Matrix{ComplexF64}
    N_k::Int
    N_b::Int
end

"""
Exact reduced map for a kernel independent of both radial momenta.

See [Microscopic kernels](@ref theory-kernels) and the exact
momentum-independent branch in [Production backend](@ref theory-production).
"""
struct MomentumIndependentProductionKernel <: AbstractProductionKernel
    block::Matrix{ComplexF64}
    N_k::Int
    N_b::Int
end

"""
Controlled truncated-SVD representation ``K \\approx U_r S_r V_r^\\dagger``
of the flattened six-index scattering map.  `relative_frobenius_residual` is
the exact tail norm computed from the discarded singular values of that
discrete matrix; it is not a continuum or physical-model error bound.

See [Controlled numerical approximations](@ref optimization-controlled).
"""
struct LowRankProductionKernel <: AbstractProductionKernel
    left::Matrix{ComplexF64}
    right_adjoint::Matrix{ComplexF64}
    N_k::Int
    N_b::Int
    rank::Int
    relative_frobenius_residual::Float64
end

"""
Static shift plans and exact kernel operators shared by all iterations.

See [Production backend](@ref theory-production),
[Field-periodic closure](@ref theory-periodicity), and
[Microscopic kernels](@ref theory-kernels).
"""
struct ProductionCache{P}
    W̃₊ᴱᵖ::EnergyShiftPlan
    W̃₋ᴱᵖ::EnergyShiftPlan
    W̃₊ᴸᴼ::EnergyShiftPlan
    W̃₋ᴸᴼ::EnergyShiftPlan
    kernels::Dict{Symbol,AbstractProductionKernel}
    fft_hilbert_plan::P
    product_hilbert_operator::Union{Nothing,Matrix{Float64}}
    estimate::ProductionMemoryEstimate
    source_grids::ModelGrids
    source_kernel_arrays::Dict{Symbol,Array{ComplexF64,6}}
    algorithms::AlgorithmOptions
end

"""
One retained scalar record from a production sweep point. `status` describes
the outer solve, while `scba_quality` preserves the independent strict,
diagnostic-approximate, unresolved, or invalid SCBA classification.

See [Production sweep, warm start, and recovery](@ref theory-production).
"""
struct ProductionSweepRecord
    temperature_K::Float64
    voltage_per_period_V::Float64
    field_V_per_m::Float64
    current_A_per_m2::Float64
    converged::Bool
    status::Symbol
    scba_quality::Symbol
    outer_iterations::Int
    final_scba_iterations::Int
    estimated_peak_bytes::Int
    wall_seconds::Float64
    metrics::Dict{Symbol,Float64}
    checkpoint::String
    warnings::Vector{Dict{String,Any}}

    function ProductionSweepRecord(
        temperature_K::Float64,
        voltage_per_period_V::Float64,
        field_V_per_m::Float64,
        current_A_per_m2::Float64,
        converged::Bool,
        status::Symbol,
        scba_quality::Symbol,
        outer_iterations::Int,
        final_scba_iterations::Int,
        estimated_peak_bytes::Int,
        wall_seconds::Float64,
        metrics::Dict{Symbol,Float64},
        checkpoint::String,
        warnings::Vector{Dict{String,Any}} = Dict{String,Any}[],
    )
        _check_scba_result_classification(converged, status, scba_quality)
        return new(
            temperature_K,
            voltage_per_period_V,
            field_V_per_m,
            current_A_per_m2,
            converged,
            status,
            scba_quality,
            outer_iterations,
            final_scba_iterations,
            estimated_peak_bytes,
            wall_seconds,
            metrics,
            checkpoint,
            warnings,
        )
    end
end

function _solution_report_metrics(solution::NEGFSolution)
    values = Dict{Symbol,Float64}()
    if !isempty(solution.outer_history)
        outer = last(solution.outer_history)
        values[:current_continuity] = outer.r_J
        values[:charge_neutrality] = outer.r_neutral
        values[:poisson] = outer.r_P
        values[:hartree] = outer.r_U
        values[:density] = outer.r_n
    end
    if !isempty(solution.scba.history)
        inner = last(solution.scba.history)
        values[:dyson] = inner.r_D
        values[:spectral_identity] = inner.r_A
        values[:keldysh] = inner.r_K
        values[:self_energy] = inner.r_Σ
        values[:normalization] = inner.r_λ
        values[:positivity] = inner.r_PSD
        values[:causality] = inner.r_caus
        values[:fft_roundoff] = inner.r_roundoff
    end
    # The final independent acceptance uses the stored, aligned state, which
    # can differ from the last outer iterate. Report its values whenever
    # available and keep the spectral sum rule visible in ordinary summaries.
    for (label, metric) in (
        (:current_continuity, :r_J),
        (:charge_neutrality, :r_neutral),
        (:poisson, :r_P),
        (:hartree, :r_U),
        (:density, :r_n),
        (:dyson, :r_D),
        (:spectral_identity, :r_A),
        (:keldysh, :r_K),
        (:self_energy, :r_Σ),
        (:normalization, :r_λ),
        (:positivity, :r_PSD),
        (:causality, :r_caus),
        (:fft_roundoff, :r_roundoff),
        (:spectral_sum, :r_sum),
        (:power_balance, :r_power),
    )
        haskey(solution.report.metrics, metric) &&
            (values[label]=solution.report.metrics[metric])
    end

    return values
end

"""
Memory-safe sweep result: summaries only, never all full NEGF states.

See [Production sweep, warm start, and recovery](@ref theory-production) and
the [state persistence contract](https://github.com/AfonenkoA/QCLNEGFRunner.jl/blob/main/docs/src/user/results.md).
"""
struct ProductionSweepResult
    records::Vector{ProductionSweepRecord}
    summary_path::String
end

function _check_production_options(options::ProductionOptions)
    _check_algorithm_options(options.algorithms)
    options.memory_budget_bytes > 0 ||
        throw(ArgumentError("memory_budget_bytes must be positive"))
    options.energy_chunk > 0 || throw(ArgumentError("energy_chunk must be positive"))
    options.hilbert_columns > 0 || throw(ArgumentError("hilbert_columns must be positive"))
    options.parallel_backend in (:threads, :blas) ||
        throw(ArgumentError("parallel_backend must be :threads or :blas"))
    options.worker_count ≥ 0 || throw(ArgumentError("worker_count cannot be negative"))
    options.residual_chunk > 0 || throw(ArgumentError("residual_chunk must be positive"))
    options.checkpoint_every_outer ≥ 0 ||
        throw(ArgumentError("checkpoint_every_outer cannot be negative"))
    options.checkpoint_every_scba ≥ 0 ||
        throw(ArgumentError("checkpoint_every_scba cannot be negative"))
    options.progress_every_scba ≥ 0 ||
        throw(ArgumentError("progress_every_scba cannot be negative"))
    options.progress_every_outer ≥ 0 ||
        throw(ArgumentError("progress_every_outer cannot be negative"))
    return options
end

"""Ask the application at a safe boundary, before allocating a snapshot.

The request contains counters only; persistence and control paths stay outside
numerics. An explicit request owns cadence. Without it the numerical cadence is
used, preserving the ordinary low-level solver contract.
"""
function _production_checkpoint_requested(
    options::ProductionOptions,
    stage::Symbol,
    iteration::Int,
)
    request = options.checkpoint_request
    if request !== nothing
        return request((; stage, iteration, outer_iteration = options.outer_iteration)) ===
               true
    end
    cadence =
        stage === :scba ? options.checkpoint_every_scba : options.checkpoint_every_outer
    return cadence > 0 && iteration % cadence == 0
end
