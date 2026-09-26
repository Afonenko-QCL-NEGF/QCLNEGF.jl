"""Numerical construction backend. Implementations consume validated physical
models but are forbidden from changing the enabled physics."""
abstract type AbstractScatteringKernelBackend end

"""E0 literal six-index kernel construction used as the test oracle."""
struct LiteralScatteringKernelBackend <: AbstractScatteringKernelBackend end

"""E1 exact production builder with threaded finite-angle blocks."""
struct ExactParallelScatteringKernelBackend <: AbstractScatteringKernelBackend end

"""Continuum angular quadrature for the impurity channel; tolerance is measured."""
struct AdaptiveDirectScatteringKernelBackend <: AbstractScatteringKernelBackend
    angular_tolerance::Float64
    function AdaptiveDirectScatteringKernelBackend(tolerance::Real = 1e-5)
        isfinite(tolerance) && 0 < tolerance < 1 ||
            throw(ArgumentError("angular tolerance must lie in (0,1)"))
        new(Float64(tolerance))
    end
end

"""E2 adaptive lookup construction with measured interpolation residuals."""
struct TabulatedScatteringKernelBackend <: AbstractScatteringKernelBackend
    options::ProductionKernelOptions
    function TabulatedScatteringKernelBackend(options::ProductionKernelOptions)
        try
            _check_production_kernel_options(options)
        catch error
            _scattering_error(:backend, "tabulated", sprint(showerror, error))
        end
        new(options)
    end
end

TabulatedScatteringKernelBackend() =
    TabulatedScatteringKernelBackend(ProductionKernelOptions())

"""Declare the E0--E2 evidence class of a numerical kernel backend."""
scattering_backend_evidence_class(::LiteralScatteringKernelBackend) = :E0
scattering_backend_evidence_class(::ExactParallelScatteringKernelBackend) = :E1
scattering_backend_evidence_class(::AdaptiveDirectScatteringKernelBackend) = :E2
scattering_backend_evidence_class(::TabulatedScatteringKernelBackend) = :E2

function scattering_backend_evidence_class(backend::AbstractScatteringKernelBackend)
    _scattering_error(
        :backend,
        "evidence_class",
        "backend $(typeof(backend)) must declare an E0, E1, or E2 class",
    )
end

"""
Numerical policy for scattering construction and contraction. `contraction`
is `:literal`, `:dense_blas`, or `:low_rank`. A nonzero low-rank tolerance or
finite rank cap is E2; exact algebraic reordering is E1.
"""
struct ScatteringNumericalPlan{B<:AbstractScatteringKernelBackend}
    kernel_backend::B
    contraction::Symbol
    low_rank_relative_tolerance::Float64
    low_rank_maximum_rank::Int
    function ScatteringNumericalPlan(
        kernel_backend::B;
        contraction::Symbol = :literal,
        low_rank_relative_tolerance::Real = 0.0,
        low_rank_maximum_rank::Integer = 0,
    ) where {B<:AbstractScatteringKernelBackend}
        contraction in (:literal, :dense_blas, :low_rank) || _scattering_error(
            :backend,
            "contraction",
            "must be literal, dense_blas, or low_rank",
        )
        0 <= low_rank_relative_tolerance < 1 ||
            _scattering_error(:backend, "low_rank_relative_tolerance", "must lie in [0,1)")
        low_rank_maximum_rank >= 0 ||
            _scattering_error(:backend, "low_rank_maximum_rank", "must be nonnegative")
        contraction === :low_rank ||
            (iszero(low_rank_relative_tolerance) && iszero(low_rank_maximum_rank)) ||
            _scattering_error(
                :backend,
                "low_rank",
                "low-rank controls require contraction=low_rank",
            )
        new{B}(
            kernel_backend,
            contraction,
            Float64(low_rank_relative_tolerance),
            Int(low_rank_maximum_rank),
        )
    end
end

"""Return the strongest E0--E3 class selected by a scattering plan/options."""
function scattering_evidence_class(plan::ScatteringNumericalPlan)
    backend_class = scattering_backend_evidence_class(plan.kernel_backend)
    backend_class in (:E0, :E1, :E2) || _scattering_error(
        :backend,
        "evidence_class",
        "numerical backend must be E0, E1, or E2",
    )
    backend_class === :E2 && return :E2
    if plan.contraction === :low_rank &&
       (!iszero(plan.low_rank_relative_tolerance) || !iszero(plan.low_rank_maximum_rank))
        return :E2
    elseif backend_class === :E1 || plan.contraction !== :literal
        return :E1
    end
    return :E0
end

"""Translate the existing public algorithm configuration into the separated
scattering numerical policy without interpreting any physical parameters."""
function scattering_numerical_plan(
    options::AlgorithmOptions,
    kernel_options::ProductionKernelOptions = ProductionKernelOptions();
    solver_backend::Symbol = options.solver_backend,
)
    try
        _check_algorithm_options(options)
    catch error
        _scattering_error(:backend, "algorithms", sprint(showerror, error))
    end
    solver_backend in (:educational, :production) ||
        _scattering_error(:backend, "solver_backend", "must be educational or production")
    solver_backend === options.solver_backend || _scattering_error(
        :backend,
        "solver_backend",
        "must match algorithms.solver_backend=$(options.solver_backend)",
    )
    backend =
        options.kernel_build === :adaptive_direct ?
        AdaptiveDirectScatteringKernelBackend(kernel_options.impurity_angular_tolerance) :
        options.kernel_build === :direct ?
        (
            solver_backend === :production ? ExactParallelScatteringKernelBackend() :
            LiteralScatteringKernelBackend()
        ) : TabulatedScatteringKernelBackend(kernel_options)
    return ScatteringNumericalPlan(
        backend;
        contraction = options.contraction,
        low_rank_relative_tolerance = options.low_rank_relative_tolerance,
        low_rank_maximum_rank = options.low_rank_maximum_rank,
    )
end

function scattering_evidence_class(
    options::AlgorithmOptions,
    kernel_options::ProductionKernelOptions = ProductionKernelOptions(),
)
    if options.retarded_real_part !== :kramers_kronig ||
       options.self_energy_structure !== :full ||
       options.transverse_momentum !== :resolved
        return :E3
    end
    return scattering_evidence_class(scattering_numerical_plan(options, kernel_options))
end

"""Typed assembly boundary between physical models and discretization data.

No numerical method may construct this context from raw unchecked
dictionaries. `ScatteringOptions` remains the documented mechanism-selection
contract; the validated model collection owns the corresponding physical
parameters.
"""
struct ScatteringAssemblyContext
    physical::PhysicalParameters
    numerical::NumericalParameters
    scales::ScaleSystem
    grids::ModelGrids
    profiles::MaterialProfiles
    basis::BasisData
end

function _require_finite_array(mechanism, field, values)
    all(value -> isfinite(real(value)) && isfinite(imag(value)), values) ||
        _scattering_error(mechanism, field, "contains non-finite values")
end

"""Validate scattering grid, profile, basis shapes and finite-value formats."""
function validate_scattering_context(context::ScatteringAssemblyContext)
    n = context.numerical
    g = context.grids
    p = context.profiles
    b = context.basis
    length(g.x) == n.N_z || _scattering_error(
        :data,
        "grids.x",
        "expected $(n.N_z) spatial nodes, got $(length(g.x))",
    )
    length(g.wˣ) == n.N_z ||
        _scattering_error(:data, "grids.wˣ", "must have the same length as grids.x")
    length(g.ε) == n.N_E ||
        _scattering_error(:data, "grids.ε", "expected $(n.N_E) energy nodes")
    length(g.wᴱ) == n.N_E ||
        _scattering_error(:data, "grids.wᴱ", "must have the same length as grids.ε")
    length(g.κ) == n.N_k ||
        _scattering_error(:data, "grids.κ", "expected $(n.N_k) radial-momentum nodes")
    length(g.wᵏ) == n.N_k ||
        _scattering_error(:data, "grids.wᵏ", "must have the same length as grids.κ")
    length(g.φ) == n.N_φ ||
        _scattering_error(:data, "grids.φ", "expected $(n.N_φ) angular nodes")
    length(g.qᶻ) == n.N_qz || _scattering_error(
        :data,
        "grids.qᶻ",
        "expected $(n.N_qz) longitudinal-momentum nodes",
    )
    length(g.wᑫᶻ) == n.N_qz ||
        _scattering_error(:data, "grids.wᑫᶻ", "must have the same length as grids.qᶻ")
    for field in fieldnames(MaterialProfiles)
        length(getfield(p, field)) == n.N_z ||
            _scattering_error(:data, "profiles.$field", "expected $(n.N_z) entries")
    end
    size(b.χ) == (n.N_z, n.N_b) || _scattering_error(
        :data,
        "basis.χ",
        "expected shape ($(n.N_z), $(n.N_b)), got $(size(b.χ))",
    )
    size(b.H₀) == (n.N_b, n.N_b) || _scattering_error(
        :data,
        "basis.H₀",
        "expected square basis matrix of order $(n.N_b)",
    )
    for (field, values) in (
        ("grids.x", g.x),
        ("grids.wˣ", g.wˣ),
        ("grids.ε", g.ε),
        ("grids.wᴱ", g.wᴱ),
        ("grids.κ", g.κ),
        ("grids.wᵏ", g.wᵏ),
        ("grids.qᶻ", g.qᶻ),
        ("grids.wᑫᶻ", g.wᑫᶻ),
        ("grids.φ", g.φ),
        ("basis.χ", b.χ),
    )
        _require_finite_array(:data, field, values)
    end
    issorted(g.x) || _scattering_error(:data, "grids.x", "spatial nodes must be sorted")
    issorted(g.ε) || _scattering_error(:data, "grids.ε", "energy nodes must be sorted")
    issorted(g.κ) ||
        _scattering_error(:data, "grids.κ", "radial-momentum nodes must be sorted")
    issorted(g.qᶻ) ||
        _scattering_error(:data, "grids.qᶻ", "longitudinal-momentum nodes must be sorted")
    all(>(0), g.wˣ) ||
        _scattering_error(:data, "grids.wˣ", "quadrature weights must be positive")
    all(>(0), g.wᴱ) ||
        _scattering_error(:data, "grids.wᴱ", "quadrature weights must be positive")
    all(>=(0), g.wᵏ) && any(>(0), g.wᵏ) || _scattering_error(
        :data,
        "grids.wᵏ",
        "radial weights must be nonnegative and not all zero",
    )
    all(>(0), g.wᑫᶻ) ||
        _scattering_error(:data, "grids.wᑫᶻ", "quadrature weights must be positive")
    return context
end



"""Validate the common six-axis normalized `KernelSet` data contract."""
function validate_scattering_kernel_set(
    kernels::KernelSet,
    numerical::NumericalParameters,
    models::AbstractVector{<:AbstractScatteringPhysicalModel};
    normalization_tolerance::Real = 64eps(Float64),
)
    isfinite(normalization_tolerance) && normalization_tolerance >= 0 || _scattering_error(
        :data,
        "normalization_tolerance",
        "must be finite and nonnegative",
    )
    expected = Set(scattering_id.(models))
    enabled = Set(kernels.enabled)
    enabled == expected || _scattering_error(
        :data,
        "kernels.enabled",
        "expected $(sort!(collect(expected))), got $(sort!(collect(enabled)))",
    )
    Set(keys(kernels.K)) == enabled || _scattering_error(
        :data,
        "kernels.K",
        "kernel keys must equal the enabled mechanism set",
    )
    Set(keys(kernels.qᴷ)) == enabled || _scattering_error(
        :data,
        "kernels.qᴷ",
        "normalization keys must equal the enabled mechanism set",
    )
    allunique(kernels.enabled) ||
        _scattering_error(:data, "kernels.enabled", "mechanism identifiers must be unique")
    expected_shape = (
        numerical.N_k,
        numerical.N_k,
        numerical.N_b,
        numerical.N_b,
        numerical.N_b,
        numerical.N_b,
    )
    for mechanism in kernels.enabled
        values = kernels.K[mechanism]
        size(values) == expected_shape || _scattering_error(
            mechanism,
            "kernel",
            "expected shape $expected_shape, got $(size(values))",
        )
        eltype(values) === ComplexF64 ||
            _scattering_error(mechanism, "kernel", "must use ComplexF64 storage")
        _require_finite_array(mechanism, "kernel", values)
        scale = kernels.qᴷ[mechanism]
        isfinite(scale) && scale > 0 ||
            _scattering_error(mechanism, "normalization", "must be finite and positive")
        isapprox(
            maximum(abs, values),
            1.0;
            rtol = normalization_tolerance,
            atol = normalization_tolerance,
        ) || _scattering_error(
            mechanism,
            "kernel",
            "normalized kernel maximum must equal one",
        )
    end
    size(kernels.Fᴸᴼ) == (numerical.N_qz, numerical.N_b, numerical.N_b) ||
        _scattering_error(
            :LO,
            "longitudinal_form_factor",
            "expected shape $((numerical.N_qz, numerical.N_b, numerical.N_b))",
        )
    _require_finite_array(:LO, "longitudinal_form_factor", kernels.Fᴸᴼ)
    return kernels
end

"""Validated kernel set, diagnostics, mechanism IDs, and evidence class."""
struct ScatteringAssemblyResult{D}
    kernels::KernelSet
    diagnostics::D
    mechanisms::Vector{Symbol}
    evidence_class::Symbol
end

"""Typed result contract returned by every scattering kernel backend."""
struct ScatteringKernelBuild{D}
    kernels::KernelSet
    diagnostics::D
end

"""Numerical-backend interface. Third-party backends implement this method
and `scattering_backend_evidence_class`; the application orchestrator never
branches on concrete backend types."""
function build_scattering_kernel_set(
    backend::AbstractScatteringKernelBackend,
    context::ScatteringAssemblyContext,
    options::ScatteringOptions,
)
    _scattering_error(
        :backend,
        "kernel_backend",
        "no kernel builder is registered for $(typeof(backend))",
    )
end

function build_scattering_kernel_set(
    ::LiteralScatteringKernelBackend,
    context::ScatteringAssemblyContext,
    options::ScatteringOptions,
)
    kernels = build_kernels(
        context.physical,
        context.numerical,
        options,
        context.scales,
        context.grids,
        context.profiles,
        context.basis,
    )
    return ScatteringKernelBuild(kernels, nothing)
end

function build_scattering_kernel_set(
    ::ExactParallelScatteringKernelBackend,
    context::ScatteringAssemblyContext,
    options::ScatteringOptions,
)
    kernels, diagnostics = build_kernels_exact_parallel(
        context.physical,
        context.numerical,
        options,
        context.scales,
        context.grids,
        context.profiles,
        context.basis,
    )
    return ScatteringKernelBuild(kernels, diagnostics)
end

function build_scattering_kernel_set(
    backend::AdaptiveDirectScatteringKernelBackend,
    context::ScatteringAssemblyContext,
    options::ScatteringOptions,
)
    kernels, diagnostics = build_kernels_exact_parallel(
        context.physical,
        context.numerical,
        options,
        context.scales,
        context.grids,
        context.profiles,
        context.basis;
        impurity_angular = :adaptive,
        angular_tolerance = backend.angular_tolerance,
    )
    return ScatteringKernelBuild(kernels, diagnostics)
end

function build_scattering_kernel_set(
    backend::TabulatedScatteringKernelBackend,
    context::ScatteringAssemblyContext,
    options::ScatteringOptions,
)
    kernels, diagnostics = build_kernels_production(
        context.physical,
        context.numerical,
        options,
        context.scales,
        context.grids,
        context.profiles,
        context.basis;
        options = backend.options,
    )
    return ScatteringKernelBuild(kernels, diagnostics)
end

"""Assemble kernels through the selected numerical backend interface."""
function assemble_scattering_kernels(
    context::ScatteringAssemblyContext,
    models::AbstractVector{<:AbstractScatteringPhysicalModel},
    plan::ScatteringNumericalPlan,
)
    validate_scattering_context(context)
    foreach(model -> validate_scattering_model(model, context.physical), models)
    allunique(scattering_id.(models)) ||
        _scattering_error(:set, "mechanisms", "mechanism identifiers must be unique")
    _models_match_physical(models, context.physical) || _scattering_error(
        :set,
        "physical_parameters",
        "mechanism values differ from the physical context",
    )
    options = scattering_options(models)
    built = build_scattering_kernel_set(plan.kernel_backend, context, options)
    built isa ScatteringKernelBuild || _scattering_error(
        :backend,
        "result",
        "backend must return ScatteringKernelBuild, got " * string(typeof(built)),
    )
    kernels = built.kernels
    diagnostics = built.diagnostics
    validate_scattering_kernel_set(kernels, context.numerical, models)
    return ScatteringAssemblyResult(
        kernels,
        diagnostics,
        scattering_id.(models),
        scattering_evidence_class(plan),
    )
end

"""
    build_configured_scattering_problem(; ...)

Composition boundary used by application-configured runs. It builds common grids and
basis exactly once, delegates scattering construction through
[`build_scattering_kernel_set`](@ref), and returns one homogeneous
`NEGFProblem` independently of the selected numerical backend.
"""
function build_configured_scattering_problem(;
    physical::PhysicalParameters,
    numerical::NumericalParameters,
    scales::ScaleSystem,
    scattering::ScatteringOptions,
    algorithms::AlgorithmOptions,
    kernel_options::ProductionKernelOptions,
    solver_backend::Symbol = algorithms.solver_backend,
    validate_static::Bool = true,
    physical_models::PhysicalModelOptions = PhysicalModelOptions(),
)
    physical = resolve_physical_models(physical, physical_models)
    grids = build_grids(physical, numerical, scales)
    profiles = build_profiles(physical, numerical, scales, grids)
    basis = build_basis(
        physical,
        numerical,
        scales,
        grids,
        profiles;
        localization = algorithms.localization,
    )
    models = scattering_models(physical, scattering)
    plan = scattering_numerical_plan(algorithms, kernel_options; solver_backend)
    context = ScatteringAssemblyContext(physical, numerical, scales, grids, profiles, basis)
    assembled = assemble_scattering_kernels(context, models, plan)
    sp = _scaled_physics(physical, scales)
    if algorithms.energy_shift === :conservative_pair
        field_shift = build_conservative_shift_pair(grids.ε, grids.wᴱ, sp.Eᵖ)
        phonon_shift = build_conservative_shift_pair(grids.ε, grids.wᴱ, sp.ħωᴸᴼ)
        W₊ᴱᵖ, W₋ᴱᵖ = field_shift.plus, field_shift.minus
        W₊ᴸᴼ, W₋ᴸᴼ = phonon_shift.plus, phonon_shift.minus
    else
        shift_builder =
            solver_backend === :production && algorithms.energy_shift === :sparse_plan ?
            build_shift_plan : build_shift_matrix
        W₊ᴱᵖ = shift_builder(grids.ε, +sp.Eᵖ)
        W₋ᴱᵖ = shift_builder(grids.ε, -sp.Eᵖ)
        W₊ᴸᴼ = shift_builder(grids.ε, +sp.ħωᴸᴼ)
        W₋ᴸᴼ = shift_builder(grids.ε, -sp.ħωᴸᴼ)
    end
    problem = NEGFProblem(
        physical,
        numerical,
        scattering,
        scales,
        grids,
        profiles,
        basis,
        assembled.kernels,
        W₊ᴱᵖ,
        W₋ᴱᵖ,
        W₊ᴸᴼ,
        W₋ᴸᴼ,
        algorithms.energy_shift === :conservative_pair ?
        :finite_volume_piecewise_constant : :nodal_linear,
    )
    augmented_kernels = add_sppa_kernel(
        problem.physical,
        problem.scales,
        problem.grids,
        problem.basis,
        problem.kernels,
        physical_models.electron_electron,
    )
    problem = NEGFProblem(
        problem.physical,
        problem.numerical,
        problem.scattering,
        problem.scales,
        problem.grids,
        problem.profiles,
        problem.basis,
        augmented_kernels,
        problem.W₊ᴱᵖ,
        problem.W₋ᴱᵖ,
        problem.W₊ᴸᴼ,
        problem.W₋ᴸᴼ,
        problem.energy_shift_discretization,
        physical_models,
    )

    if validate_static
        report = validate_problem(problem)
        report.passed || throw(
            ArgumentError(
                "static configured model validation failed: " * join(report.messages, "; "),
            ),
        )
    end
    return (
        problem = problem,
        kernel_diagnostics = assembled.diagnostics,
        scattering_evidence_class = assembled.evidence_class,
    )
end
