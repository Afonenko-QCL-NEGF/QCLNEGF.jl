"""
    ProductionKernelOptions

Controls the *construction* of production scattering tensors.  LO-phonon and
ionized-impurity tensors are tabulated as functions of the transferred radial
wave number and interpolated linearly.  The table is refined until direct,
off-grid evaluations satisfy `relative_tolerance`, or until
`maximum_lookup_nodes` is reached.

`relative_tolerance` bounds all three direct diagnostics: globally scaled
entrywise error, worst pointwise block-norm error, and selected angular-average
error.  These are measured interpolation residuals, not a certified continuum
error bound.  `strict=true` refuses a table if any diagnostic misses the
requested residual.  IFR, acoustic, and alloy kernels use exact separable
constructions and do not use this tolerance.

See [Microscopic kernels](@ref theory-kernels) and
[controlled production kernel construction](@ref theory-production).
"""
Base.@kwdef struct ProductionKernelOptions
    initial_lookup_nodes::Int = 65
    maximum_lookup_nodes::Int = 513
    relative_tolerance::Float64 = 2e-4
    lookup_power::Float64 = 2.0
    validation_fractions::NTuple{3,Float64} = (0.25, 0.50, 0.75)
    angular_validation_pairs::Int = 3
    impurity_angular_tolerance::Float64 = 1e-5
    strict::Bool = true
end

"""
Measured convergence record for one microscopic kernel family.

See [Microscopic kernels](@ref theory-kernels) and
[production kernel diagnostics](@ref theory-production).
"""
struct KernelInterpolationDiagnostic
    mechanism::Symbol
    construction::Symbol
    lookup_nodes::Int
    direct_q_evaluations::Int
    sampled_global_relative_residual::Float64
    sampled_pointwise_relative_residual::Float64
    sampled_angular_relative_residual::Float64
    node_history::Vector{Int}
    residual_history::Vector{Float64}
    accepted::Bool
    angular_quadrature_checked::Bool
    angular_quadrature_relative_error::Union{Nothing,Float64}
end

# Older constructors only measured table/algebra residuals. In particular a
# zero interpolation residual never certifies angular continuum convergence.
KernelInterpolationDiagnostic(
    mechanism,
    construction,
    lookup_nodes,
    direct_q_evaluations,
    global_residual,
    pointwise_residual,
    angular_residual,
    node_history,
    residual_history,
    accepted,
) = KernelInterpolationDiagnostic(
    mechanism,
    construction,
    lookup_nodes,
    direct_q_evaluations,
    global_residual,
    pointwise_residual,
    angular_residual,
    node_history,
    residual_history,
    accepted,
    false,
    nothing,
)

"""
Diagnostics returned together with a production [`KernelSet`](@ref).

See [Microscopic kernels](@ref theory-kernels) and
[production kernel diagnostics](@ref theory-production).
"""
struct ProductionKernelDiagnostics
    mechanisms::Dict{Symbol,KernelInterpolationDiagnostic}
    worst_measured_relative_residual::Float64
    all_accepted::Bool
end

function Base.show(io::IO, diagnostic::KernelInterpolationDiagnostic)
    print(
        io,
        "KernelInterpolationDiagnostic(",
        diagnostic.mechanism,
        ", construction=",
        diagnostic.construction,
        ", nodes=",
        diagnostic.lookup_nodes,
        ", measured_residual=",
        round(
            max(
                diagnostic.sampled_global_relative_residual,
                diagnostic.sampled_pointwise_relative_residual,
                diagnostic.sampled_angular_relative_residual,
            );
            sigdigits = 4,
        ),
        ", accepted=",
        diagnostic.accepted,
        ")",
    )
end

function _check_production_kernel_options(options::ProductionKernelOptions)
    options.initial_lookup_nodes >= 3 ||
        throw(ArgumentError("initial_lookup_nodes must be at least three"))
    options.maximum_lookup_nodes >= options.initial_lookup_nodes || throw(
        ArgumentError("maximum_lookup_nodes must not be smaller than the initial table"),
    )
    0 < options.relative_tolerance < 1 ||
        throw(ArgumentError("relative_tolerance must lie in (0,1)"))
    options.lookup_power >= 1 || throw(ArgumentError("lookup_power must be at least one"))
    all(f -> 0 < f < 1, options.validation_fractions) ||
        throw(ArgumentError("validation fractions must lie strictly inside (0,1)"))
    options.angular_validation_pairs >= 0 ||
        throw(ArgumentError("angular_validation_pairs cannot be negative"))
    isfinite(options.impurity_angular_tolerance) &&
    0 < options.impurity_angular_tolerance < 1 ||
        throw(ArgumentError("impurity_angular_tolerance must lie in (0,1)"))
    return options
end

_kernel_pair_index(a::Int, c::Int, Nb::Int) = a + Nb * (c - 1)

function _pair_covariance_to_block(pair::AbstractMatrix{<:Number}, Nb::Int)
    size(pair) == (Nb^2, Nb^2) ||
        throw(DimensionMismatch("pair covariance has the wrong size"))
    # `pair[(a,c),(b,d)]` -> stored four-index order `(a,c,d,b)`.
    return ComplexF64.(permutedims(reshape(pair, Nb, Nb, Nb, Nb), (1, 2, 4, 3)))
end

"""
Create an optimized direct evaluator `K(qbar)` for the donor covariance.
Only grid cells with nonzero donor weight are retained; this is algebraically
exact.  The remaining contractions are dense matrix products over the
longitudinal quadrature and donor coordinates.
"""
function _production_impurity_evaluator(
    p::PhysicalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    profiles::MaterialProfiles,
    basis::BasisData,
)
    Nb = size(basis.χ, 2)
    Nx = length(grids.x)
    Npair = Nb^2
    donor_weights_all = grids.wˣ .* profiles.Nᴰ ./ s.L₀_m^2
    donors = findall(!iszero, donor_weights_all)
    isempty(donors) && throw(ArgumentError("impurity kernel has no ionized donors"))
    donor_weights = donor_weights_all[donors]

    pair_density = Matrix{ComplexF64}(undef, Nx, Npair)
    for c = 1:Nb, a = 1:Nb, g = 1:Nx
        pair_density[g, _kernel_pair_index(a, c, Nb)] =
            grids.wˣ[g] * conj(basis.χ[g, a]) * basis.χ[g, c]
    end
    distances = [abs(grids.x[g] - grids.x[gD]) for g = 1:Nx, gD in donors]
    eC = _electron_charge()
    ε0 = _epsilon_zero()
    qs = _inverse_metres(p.q_s)
    L0 = s.L₀_m

    function evaluate(qbar::Real)
        qbar >= 0 || throw(DomainError(qbar, "transferred wave number is negative"))
        q = Float64(qbar) / L0
        prefactor = eC^2 / (2ε0 * p.ε_s * (q + qs))
        decay = exp.(-Float64(qbar) .* distances)
        # V[(a,c),D] has exactly the same cell-centred quadrature as the
        # educational four-loop implementation.
        V = prefactor .* (transpose(pair_density) * decay)
        weighted_V = V .* reshape(donor_weights, 1, :)
        covariance = weighted_V * adjoint(V)
        return _pair_covariance_to_block(covariance, Nb) ./ _K₀(s)
    end
    return evaluate
end

"""
Create an optimized direct evaluator for the screened bulk Fröhlich tensor.
The `q_z` trapezoid is unchanged; only its four basis indices are flattened
into a covariance matrix before the BLAS product.
"""
function _production_lo_evaluator(
    p::PhysicalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    basis::BasisData,
)
    Nb = size(basis.χ, 2)
    Nq = length(grids.qᶻ)
    Fᴸᴼ = zeros(ComplexF64, Nq, Nb, Nb)
    for h = 1:Nq, a = 1:Nb, c = 1:Nb, g in eachindex(grids.x)
        Fᴸᴼ[h, a, c] +=
            grids.wˣ[g] *
            conj(basis.χ[g, a]) *
            exp(im * grids.qᶻ[h] * grids.x[g]) *
            basis.χ[g, c]
    end
    Fpair = reshape(Fᴸᴼ, Nq, Nb^2)
    eC = _electron_charge()
    ε0 = _epsilon_zero()
    ħω = Float64(ustrip(u"J", uconvert(u"J", p.ħωᴸᴼ)))
    prefactor = eC^2 * ħω / (2ε0) * (1 / p.ε_∞ - 1 / p.ε_s)
    qscreen = scaled_value(p.qᴸᴼ_s, s, :wavenumber)

    function evaluate(qbar::Real)
        qbar >= 0 || throw(DomainError(qbar, "transferred wave number is negative"))
        factors = [
            grids.wᑫᶻ[h] * s.L₀_m / (2π * (Float64(qbar)^2 + grids.qᶻ[h]^2 + qscreen^2)) for h = 1:Nq
        ]
        # No conjugation on the first vertex: the literal tensor contains
        # F[a,c] * conj(F[b,d]).
        covariance = transpose(Fpair) * (factors .* conj.(Fpair))
        return prefactor .* _pair_covariance_to_block(covariance, Nb) ./ _K₀(s)
    end
    return evaluate, Fᴸᴼ
end

function _lookup_block(nodes::AbstractVector{<:Real}, values::AbstractVector, q::Real)
    qlo, qhi = first(nodes), last(nodes)
    tolerance = 64eps(Float64) * max(abs(qlo), abs(qhi), 1.0)
    qlo - tolerance <= q <= qhi + tolerance ||
        throw(DomainError(q, "lookup query lies outside the tabulated interval"))
    qclamped = clamp(Float64(q), Float64(qlo), Float64(qhi))
    if qclamped == qhi
        return values[end]
    end
    j = clamp(searchsortedlast(nodes, qclamped), 1, length(nodes) - 1)
    θ = (qclamped - nodes[j]) / (nodes[j+1] - nodes[j])
    return (1 - θ) .* values[j] .+ θ .* values[j+1]
end

function _cached_kernel_value!(cache::Dict{Float64,Array{ComplexF64,4}}, evaluate, q::Real)
    qkey = Float64(q)
    return get!(cache, qkey) do
        ComplexF64.(evaluate(qkey))
    end
end

function _validation_metrics(nodes, values, evaluate, cache, fractions::Tuple)
    exact_values = Array{ComplexF64,4}[]
    approximate_values = Array{ComplexF64,4}[]
    for j = 1:(length(nodes)-1), fraction in fractions
        q = nodes[j] + fraction * (nodes[j+1] - nodes[j])
        push!(exact_values, _cached_kernel_value!(cache, evaluate, q))
        push!(approximate_values, _lookup_block(nodes, values, q))
    end
    global_scale = maximum(maximum(abs, value) for value in (values..., exact_values...))
    global_scale > 0 || return 0.0, 0.0
    global_residual =
        maximum(
            maximum(abs, exact - approximate) for
            (exact, approximate) in zip(exact_values, approximate_values)
        ) / global_scale
    pointwise_residual = maximum(
        norm(exact - approximate) /
        max(norm(exact), sqrt(length(exact)) * global_scale * 1e-12) for
        (exact, approximate) in zip(exact_values, approximate_values)
    )
    return global_residual, pointwise_residual
end

function _selected_momentum_pairs(Nk::Int, requested::Int)
    requested == 0 && return Tuple{Int,Int}[]
    middle = cld(Nk, 2)
    candidates = [
        (1, Nk),
        (middle, middle),
        (Nk, Nk),
        (1, 1),
        (1, middle),
        (middle, Nk),
        (Nk, 1),
        (middle, 1),
        (Nk, middle),
    ]
    unique!(candidates)
    return candidates[1:min(requested, length(candidates))]
end

function _angular_validation_metrics(
    nodes,
    values,
    evaluate,
    cache,
    grids::ModelGrids,
    requested_pairs::Int,
)
    pairs = _selected_momentum_pairs(length(grids.κ), requested_pairs)
    isempty(pairs) && return 0.0
    exact_blocks = Array{ComplexF64,4}[]
    approximate_blocks = Array{ComplexF64,4}[]
    for (m, mp) in pairs
        exact = zeros(ComplexF64, size(values[1]))
        approximate = similar(exact)
        fill!(approximate, 0)
        for φ in grids.φ
            q = sqrt(
                max(grids.κ[m]^2 + grids.κ[mp]^2 - 2grids.κ[m] * grids.κ[mp] * cos(φ), 0.0),
            )
            exact .+= _cached_kernel_value!(cache, evaluate, q)
            approximate .+= _lookup_block(nodes, values, q)
        end
        exact ./= length(grids.φ)
        approximate ./= length(grids.φ)
        push!(exact_blocks, exact)
        push!(approximate_blocks, approximate)
    end
    scale = maximum(maximum(abs, block) for block in exact_blocks)
    scale > 0 || return 0.0
    return maximum(
        maximum(abs, exact - approximate) for
        (exact, approximate) in zip(exact_blocks, approximate_blocks)
    ) / scale
end

function _assemble_lookup_kernel(nodes, values, grids::ModelGrids, Nb::Int)
    Nk = length(grids.κ)
    K = zeros(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
    Base.Threads.@threads :static for linear = 1:(Nk*Nk)
        local m = (linear - 1) % Nk + 1
        local mp = (linear - 1) ÷ Nk + 1
        local destination = view(K, m, mp, :, :, :, :)
        for φ in grids.φ
            local q = sqrt(
                max(grids.κ[m]^2 + grids.κ[mp]^2 - 2grids.κ[m] * grids.κ[mp] * cos(φ), 0.0),
            )
            destination .+= _lookup_block(nodes, values, q)
        end
        destination ./= length(grids.φ)
    end
    return K
end

"""
Assemble the finite-angle tensor without interpolation. The evaluator uses the
same quadratures as the formula-level implementation; only independent
`(k,k')` blocks and basis contractions are reordered for production threads.
"""
function _assemble_exact_parallel_kernel(evaluate, grids::ModelGrids, Nb::Int)
    Nk = length(grids.κ)
    K = zeros(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
    Base.Threads.@threads :static for linear = 1:(Nk*Nk)
        local m = (linear - 1) % Nk + 1
        local mp = (linear - 1) ÷ Nk + 1
        local destination = view(K, m, mp, :, :, :, :)
        for φ in grids.φ
            local q = sqrt(
                max(grids.κ[m]^2 + grids.κ[mp]^2 - 2grids.κ[m] * grids.κ[mp] * cos(φ), 0.0),
            )
            destination .+= evaluate(q)
        end
        destination ./= length(grids.φ)
    end
    return K
end

function _tabulated_angular_kernel(
    mechanism::Symbol,
    evaluate,
    grids::ModelGrids,
    Nb::Int,
    options::ProductionKernelOptions,
)
    qmaximum = 2maximum(grids.κ)
    qmaximum > 0 || throw(ArgumentError("radial grid has zero extent"))
    cache = Dict{Float64,Array{ComplexF64,4}}()
    node_history = Int[]
    residual_history = Float64[]
    nodes_count = options.initial_lookup_nodes
    final_nodes = Float64[]
    final_values = Array{ComplexF64,4}[]
    final_global = Inf
    final_pointwise = Inf
    final_angular = Inf
    accepted = false

    while true
        # Screened Coulomb and Fröhlich kernels have their largest curvature
        # close to q=0.  Uniform auxiliary nodes followed by q=qmax*t^p keep
        # the table nested under refinement and place resolution where it is
        # needed, while interpolation itself remains linear in physical q.
        nodes = [
            qmaximum * ((j - 1) / (nodes_count - 1))^options.lookup_power for
            j = 1:nodes_count
        ]
        values = [_cached_kernel_value!(cache, evaluate, q) for q in nodes]
        global_residual, pointwise_residual = _validation_metrics(
            nodes,
            values,
            evaluate,
            cache,
            options.validation_fractions,
        )
        angular_residual = _angular_validation_metrics(
            nodes,
            values,
            evaluate,
            cache,
            grids,
            options.angular_validation_pairs,
        )
        measured = maximum((global_residual, pointwise_residual, angular_residual))
        push!(node_history, nodes_count)
        push!(residual_history, measured)
        final_nodes, final_values = nodes, values
        final_global, final_pointwise = global_residual, pointwise_residual
        final_angular = angular_residual
        accepted = isfinite(measured) && measured <= options.relative_tolerance
        (accepted || nodes_count == options.maximum_lookup_nodes) && break
        nodes_count = min(2 * (nodes_count - 1) + 1, options.maximum_lookup_nodes)
    end

    diagnostic = KernelInterpolationDiagnostic(
        mechanism,
        :piecewise_linear_q_lookup,
        length(final_nodes),
        length(cache),
        final_global,
        final_pointwise,
        final_angular,
        node_history,
        residual_history,
        accepted,
    )
    if !accepted && options.strict
        throw(
            ErrorException(
                "$mechanism q-lookup failed the measured " *
                "relative tolerance $(options.relative_tolerance); " *
                "nodes=$(length(final_nodes)), residual=$(maximum((final_global, final_pointwise, final_angular)))",
            ),
        )
    end
    K = _assemble_lookup_kernel(final_nodes, final_values, grids, Nb)
    return K, diagnostic
end

"""
Exact finite-angular-grid IFR construction.  The basis covariance is formed
once, while the Gaussian spectrum is averaged as a scalar for every `(k,k')`
pair.  This is algebraically identical to `interface_roughness_kernel` for the
same `grids.φ`; no continuum-angle or Bessel-function approximation is used.
"""
function _interface_roughness_kernel_separable(
    p::PhysicalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    basis::BasisData,
)
    Nb = size(basis.χ, 2)
    Δ = _metres(p.Δᴵᶠᴿ)
    Λ = _metres(p.Λᴵᶠᴿ)
    L0 = s.L₀_m
    covariance = zeros(ComplexF64, Nb^2, Nb^2)
    for zint in p.interfaces
        xint = _metres(zint) / L0
        χint = [_interface_value(view(basis.χ, :, a), grids.x, xint) for a = 1:Nb]
        ΔEc = Float64(ustrip(u"J", _band_jump(p, zint)))
        D = [ΔEc * conj(χint[a]) * χint[c] / L0 for a = 1:Nb, c = 1:Nb]
        vector = reshape(D, Nb^2)
        covariance .+= vector * adjoint(vector)
    end
    block = _pair_covariance_to_block(covariance, Nb) ./ _K₀(s)
    Nk = length(grids.κ)
    K = zeros(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
    ratio2 = (Λ / L0)^2
    amplitude = π * Δ^2 * Λ^2
    Base.Threads.@threads :static for linear = 1:(Nk*Nk)
        local m = (linear - 1) % Nk + 1
        local mp = (linear - 1) ÷ Nk + 1
        local angular_factor = 0.0
        for φ in grids.φ
            local qbar2 =
                max(grids.κ[m]^2 + grids.κ[mp]^2 - 2grids.κ[m] * grids.κ[mp] * cos(φ), 0.0)
            angular_factor += exp(-ratio2 * qbar2 / 4)
        end
        angular_factor *= amplitude / length(grids.φ)
        view(K, m, mp, :, :, :, :) .= angular_factor .* block
    end
    return K
end

function _exact_kernel_diagnostic(mechanism::Symbol, construction::Symbol)
    return KernelInterpolationDiagnostic(
        mechanism,
        construction,
        0,
        0,
        0.0,
        0.0,
        0.0,
        Int[],
        Float64[],
        true,
    )
end

function _normalize_production_kernels(
    raw::Dict{Symbol,Array{ComplexF64,6}},
    Fᴸᴼ::Array{ComplexF64,3},
)
    requested = [name for name in _SCATTERING_MECHANISM_ORDER if haskey(raw, name)]
    enabled = Symbol[]
    qᴷ = Dict{Symbol,Float64}()
    normalized = Dict{Symbol,Array{ComplexF64,6}}()
    for name in requested
        q = maximum(abs, raw[name])
        isfinite(q) || throw(DomainError(q, "kernel $name is non-finite"))
        if iszero(q)
            @info "zero microscopic kernel is declared disabled" mechanism=name
            continue
        end
        q > 0 || throw(DomainError(q, "kernel scale must be positive"))
        push!(enabled, name)
        qᴷ[name] = q
        # `raw` is private construction storage and is not returned or shared.
        # Normalize it in place so the static production stage never retains a
        # second complete set of six-axis tensors merely to rescale them.
        raw[name] ./= q
        normalized[name] = raw[name]
    end
    return KernelSet(normalized, qᴷ, Fᴸᴼ, enabled)
end

"""
    build_kernels_production(physical, numerical, scattering, scales,
                             grids, profiles, basis; options)

Build the same full six-axis `KernelSet` contract as [`build_kernels`](@ref),
but avoid repeated expensive longitudinal integrations:

* LO and impurity use adaptively refined, linearly interpolated `q` tables;
* IFR uses an exact separable discrete-angle construction;
* acoustic and alloy use their exact momentum-independent constructions.

Returns `(kernel_set, diagnostics)`.  The diagnostic residuals are direct
off-grid measurements for this concrete basis and grid.  They are neither an
analytic interpolation bound nor a substitute for a final observable-level
convergence run with a tighter table tolerance.

See [Microscopic kernels](@ref theory-kernels) for the physical tensors and
[Production backend](@ref theory-production) for the controlled construction
and residual interpretation.
"""
function build_kernels_production(
    p::PhysicalParameters,
    n::NumericalParameters,
    scattering::ScatteringOptions,
    s::ScaleSystem,
    grids::ModelGrids,
    profiles::MaterialProfiles,
    basis::BasisData;
    options::ProductionKernelOptions = ProductionKernelOptions(),
)
    _check_production_kernel_options(options)
    Nb = size(basis.χ, 2)
    raw = Dict{Symbol,Array{ComplexF64,6}}()
    diagnostics = Dict{Symbol,KernelInterpolationDiagnostic}()
    Fᴸᴼ = zeros(ComplexF64, length(grids.qᶻ), Nb, Nb)

    if scattering.impurity
        evaluator = _production_impurity_evaluator(p, s, grids, profiles, basis)
        raw[:impurity], diagnostics[:impurity] =
            _tabulated_angular_kernel(:impurity, evaluator, grids, Nb, options)
    end
    if scattering.IFR
        raw[:IFR] = _interface_roughness_kernel_separable(p, s, grids, basis)
        diagnostics[:IFR] = _exact_kernel_diagnostic(:IFR, :exact_separable_angle)
    end
    if scattering.acoustic
        raw[:acoustic] = acoustic_kernel(p, s, grids, basis)
        diagnostics[:acoustic] =
            _exact_kernel_diagnostic(:acoustic, :exact_momentum_independent)
    end
    if scattering.alloy
        p.ΔV_alloy === nothing &&
            throw(ArgumentError("alloy enabled but ΔV_alloy is missing"))
        p.Ω₀ === nothing && throw(ArgumentError("alloy enabled but Ω₀ is missing"))
        raw[:alloy] = alloy_kernel(p, s, grids, profiles, basis)
        diagnostics[:alloy] = _exact_kernel_diagnostic(:alloy, :exact_momentum_independent)
    end
    if scattering.LO
        evaluator, Fᴸᴼ = _production_lo_evaluator(p, s, grids, basis)
        raw[:LO], diagnostics[:LO] =
            _tabulated_angular_kernel(:LO, evaluator, grids, Nb, options)
    end

    kernel_set = _normalize_production_kernels(raw, Fᴸᴼ)
    worst =
        isempty(diagnostics) ? 0.0 :
        maximum(
            maximum((
                d.sampled_global_relative_residual,
                d.sampled_pointwise_relative_residual,
                d.sampled_angular_relative_residual,
            )) for d in values(diagnostics)
        )
    result_diagnostics = ProductionKernelDiagnostics(
        diagnostics,
        worst,
        all(d.accepted for d in values(diagnostics)),
    )
    return kernel_set, result_diagnostics
end

"""
Build exact direct scattering tensors with the optimized algebraic evaluators.
No lookup table, interpolation or model reduction is used. This is the E1
production counterpart of the unchanged literal educational builder.
"""
function build_kernels_exact_parallel(
    p::PhysicalParameters,
    n::NumericalParameters,
    scattering::ScatteringOptions,
    s::ScaleSystem,
    grids::ModelGrids,
    profiles::MaterialProfiles,
    basis::BasisData;
    impurity_angular::Symbol = :legacy_uniform,
    angular_tolerance::Real = 1e-5,
)
    impurity_angular in (:legacy_uniform, :adaptive) ||
        throw(ArgumentError("unsupported impurity angular quadrature"))
    Nb = size(basis.χ, 2)
    raw = Dict{Symbol,Array{ComplexF64,6}}()
    diagnostics = Dict{Symbol,KernelInterpolationDiagnostic}()
    Fᴸᴼ = zeros(ComplexF64, length(grids.qᶻ), Nb, Nb)

    if scattering.impurity
        evaluator = _production_impurity_evaluator(p, s, grids, profiles, basis)
        if impurity_angular === :legacy_uniform
            raw[:impurity] = _assemble_exact_parallel_kernel(evaluator, grids, Nb)
            diagnostics[:impurity] =
                _exact_kernel_diagnostic(:impurity, :exact_parallel_finite_angle)
        else
            Nk = length(grids.κ)
            tensor = zeros(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
            errors, evaluations = zeros(Float64, Nk*Nk), zeros(Int, Nk*Nk)
            accepted = fill(false, Nk*Nk)
            Base.Threads.@threads :static for linear = 1:(Nk*Nk)
                m, mp = (linear-1)%Nk+1, (linear-1)÷Nk+1
                result = adaptive_angular_average(
                    evaluator,
                    grids.κ[m],
                    grids.κ[mp];
                    relative_tolerance = angular_tolerance,
                )
                tensor[m, mp, :, :, :, :] .= result.value
                errors[linear], evaluations[linear], accepted[linear] =
                    result.estimated_relative_error, result.evaluations, result.accepted
            end
            all(accepted) || throw(
                ErrorException(
                    "adaptive impurity angular integration failed; estimated relative error=$(maximum(errors))",
                ),
            )
            raw[:impurity] = tensor
            diagnostics[:impurity] = KernelInterpolationDiagnostic(
                :impurity,
                :adaptive_transformed_simpson,
                0,
                sum(evaluations),
                0.0,
                0.0,
                maximum(errors),
                Int[],
                [maximum(errors)],
                true,
                true,
                maximum(errors),
            )
        end
    end
    if scattering.IFR
        raw[:IFR] = _interface_roughness_kernel_separable(p, s, grids, basis)
        diagnostics[:IFR] = _exact_kernel_diagnostic(:IFR, :exact_parallel_separable_angle)
    end
    if scattering.acoustic
        raw[:acoustic] = acoustic_kernel(p, s, grids, basis)
        diagnostics[:acoustic] =
            _exact_kernel_diagnostic(:acoustic, :exact_momentum_independent)
    end
    if scattering.alloy
        p.ΔV_alloy === nothing &&
            throw(ArgumentError("alloy enabled but ΔV_alloy is missing"))
        p.Ω₀ === nothing && throw(ArgumentError("alloy enabled but Ω₀ is missing"))
        raw[:alloy] = alloy_kernel(p, s, grids, profiles, basis)
        diagnostics[:alloy] = _exact_kernel_diagnostic(:alloy, :exact_momentum_independent)
    end
    if scattering.LO
        evaluator, Fᴸᴼ = _production_lo_evaluator(p, s, grids, basis)
        raw[:LO] = _assemble_exact_parallel_kernel(evaluator, grids, Nb)
        diagnostics[:LO] = _exact_kernel_diagnostic(:LO, :exact_parallel_finite_angle)
    end

    kernels = _normalize_production_kernels(raw, Fᴸᴼ)
    worst =
        isempty(diagnostics) ? 0.0 :
        maximum(d.sampled_angular_relative_residual for d in values(diagnostics))
    return kernels, ProductionKernelDiagnostics(diagnostics, worst, true)
end

function _preflight_production_memory(
    numerical::NumericalParameters,
    scattering::ScatteringOptions,
    production_options::ProductionOptions,
)
    requested = count(
        identity,
        (
            scattering.LO,
            scattering.acoustic,
            scattering.impurity,
            scattering.IFR,
            scattering.alloy,
        ),
    )
    dense = count(identity, (scattering.LO, scattering.impurity, scattering.IFR))
    lo_kernel = scattering.LO ? :dense : :absent
    estimate = estimate_production_memory(
        numerical,
        requested;
        dense_mechanism_count = dense,
        lo_kernel,
        options = production_options,
    )
    estimate.peak_bytes <= production_options.memory_budget_bytes || throw(
        ArgumentError(
            "estimated production peak $(estimate.peak_bytes) bytes exceeds " *
            "the resolved execution budget " *
            "$(production_options.memory_budget_bytes) bytes",
        ),
    )
    return estimate
end

"""
    build_problem_production(; physical=reference_parameters(),
        numerical=reference_production_numerics(), scattering=default_scattering(),
        scales=ScaleSystem(), production_options=ProductionOptions(),
        kernel_options=ProductionKernelOptions())

Production static-stage counterpart of [`build_problem`](@ref). A
conservative RAM preflight is performed before the spatial basis and kernels
are allocated. The spatial
discretization, basis, complete six-axis kernel contract, and shift-matrix
semantics are unchanged.  Only kernel construction is replaced by
[`build_kernels_production`](@ref).  Returning diagnostics alongside the
problem makes the measured interpolation residual part of every production
run rather than an unrecorded build setting.

The returned named tuple contains `problem`, `kernel_diagnostics`, and
`memory_estimate`, so the measured interpolation error cannot become an
unrecorded build setting. With `energy_shift=:sparse_plan`, `NEGFProblem`
stores the exact `O(N_E)` interpolation operators directly. The educational
and explicitly dense paths retain the dense matrix oracle.

See [Production backend](@ref theory-production),
[numerical cost and memory](@ref theory-cost), and
[the implementation plan](@ref implementation-plan).
"""
function build_problem_production(;
    physical::PhysicalParameters = reference_parameters(),
    numerical::NumericalParameters = reference_production_numerics(),
    scattering::ScatteringOptions = default_scattering(),
    scales::ScaleSystem = ScaleSystem(),
    production_options::ProductionOptions = ProductionOptions(),
    kernel_options::ProductionKernelOptions = ProductionKernelOptions(),
    localization::Symbol = :pzp,
    validate_static::Bool = true,
    physical_models::PhysicalModelOptions = PhysicalModelOptions(),
)
    physical = resolve_physical_models(physical, physical_models)
    estimate = _preflight_production_memory(numerical, scattering, production_options)
    grids = build_grids(physical, numerical, scales)
    profiles = build_profiles(physical, numerical, scales, grids)
    basis = build_basis(physical, numerical, scales, grids, profiles; localization)
    kernels, diagnostics = build_kernels_production(
        physical,
        numerical,
        scattering,
        scales,
        grids,
        profiles,
        basis;
        options = kernel_options,
    )
    sp = _scaled_physics(physical, scales)
    shift_builder =
        production_options.algorithms.energy_shift === :sparse_plan ? build_shift_plan :
        build_shift_matrix
    W₊ᴱᵖ = shift_builder(grids.ε, +sp.Eᵖ)
    W₋ᴱᵖ = shift_builder(grids.ε, -sp.Eᵖ)
    W₊ᴸᴼ = shift_builder(grids.ε, +sp.ħωᴸᴼ)
    W₋ᴸᴼ = shift_builder(grids.ε, -sp.ħωᴸᴼ)
    problem = NEGFProblem(
        physical,
        numerical,
        scattering,
        scales,
        grids,
        profiles,
        basis,
        kernels,
        W₊ᴱᵖ,
        W₋ᴱᵖ,
        W₊ᴸᴼ,
        W₋ᴸᴼ,
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
                "static production model validation failed: " * join(report.messages, "; "),
            ),
        )
    end
    return (problem = problem, kernel_diagnostics = diagnostics, memory_estimate = estimate)
end
