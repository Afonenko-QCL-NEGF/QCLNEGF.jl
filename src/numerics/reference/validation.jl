_relative_norm(A, B; floor = 1e-14) = norm(A - B) / (norm(B) + floor)

function _maximum_block_residual(X::AbstractArray{<:Number,4}, transform)
    result = 0.0
    for e in axes(X, 1), m in axes(X, 2)
        block = _matrix_block(X, e, m)
        result = max(result, transform(block))
    end
    return result
end

function _negative_part(
    X::AbstractMatrix;
    floor = 1e-14,
    construction_scale = 1.0,
    construction_error = 0.0,
)
    all(isfinite, X) || return Inf
    H = Hermitian((X + X') / 2)
    λmin = eigmin(H)
    matrix_norm = opnorm(Matrix(H))
    # Independent reference eigensolve/norm; shared contract only, no optimized kernel.
    budget =
        psd_error_budget(matrix_norm, size(X, 1); construction_scale, construction_error)
    return max(0.0, -λmin-budget) / max(matrix_norm, floatmin(Float64))
end

function _positive_imaginary_part(Σᴿ::AbstractMatrix; floor = 1e-14)
    ImΣ = Hermitian((Σᴿ - Σᴿ') / (2im))
    return max(0.0, eigmax(ImΣ)) / (opnorm(Σᴿ) + floor)
end

function _dyson_residual(ε, h, Σᴿ, Gᴿ)
    NE, Nk, Nb, _ = size(Gᴿ)
    I_b = Matrix{ComplexF64}(I, Nb, Nb)
    r = 0.0
    for e = 1:NE, m = 1:Nk
        D = ε[e] * I_b - Matrix(view(h, m, :, :)) - _matrix_block(Σᴿ, e, m)
        GR = _matrix_block(Gᴿ, e, m)
        numerator = norm(D * GR - I_b)
        denominator = norm(D) * norm(GR) + sqrt(Nb)
        r = max(r, numerator / denominator)
    end
    return r
end

function _spectral_residual(green::GreenState, total::SelfEnergyFamily)
    Γ = im .* (total.Σᵍ .- total.Σˡ)
    numerator = 0.0
    denominator = 0.0
    for e in axes(green.Gᴿ, 1), m in axes(green.Gᴿ, 2)
        GR = _matrix_block(green.Gᴿ, e, m)
        difference = _matrix_block(green.A, e, m) - GR * _matrix_block(Γ, e, m) * GR'
        numerator += norm(difference)^2
        denominator += norm(_matrix_block(green.A, e, m))^2
    end
    return sqrt(numerator) / (sqrt(denominator) + 1e-14)
end

function _keldysh_residual(green::GreenState, candidate::SelfEnergyFamily)
    numerator = 0.0
    denominator = 0.0
    for e in axes(green.Gᴿ, 1), m in axes(green.Gᴿ, 2)
        GR = _matrix_block(green.Gᴿ, e, m)
        predicted = GR * _matrix_block(candidate.Σˡ, e, m) * GR'
        difference = _matrix_block(green.Gˡ, e, m) - predicted
        numerator += norm(difference)^2
        denominator += norm(_matrix_block(green.Gˡ, e, m))^2
    end
    return sqrt(numerator) / (sqrt(denominator) + 1e-14)
end

function _selfenergy_residual(old::SelfEnergyFamily, candidate::SelfEnergyFamily)
    result = 0.0
    for name in (:Σᴿ, :Σˡ, :Σᵍ)
        Xold = getfield(old, name)
        Xnew = getfield(candidate, name)
        result = max(result, norm(Xnew - Xold) / (norm(Xnew) + 1e-14))
    end
    return result
end

function _green_positivity(green::GreenState, total::SelfEnergyFamily)
    Γ = im .* (total.Σᵍ .- total.Σˡ)
    r = 0.0
    for e in axes(green.Gᴿ, 1), m in axes(green.Gᴿ, 2)
        GR = _matrix_block(green.Gᴿ, e, m)
        r = max(
            r,
            _negative_part(_matrix_block(green.A, e, m)),
            _negative_part(-im .* _matrix_block(green.Gˡ, e, m)),
            _negative_part(im .* _matrix_block(green.Gᵍ, e, m)),
            _negative_part(-im .* (GR * _matrix_block(total.Σˡ, e, m) * GR')),
            _negative_part(im .* (GR * _matrix_block(total.Σᵍ, e, m) * GR')),
            _negative_part(_matrix_block(Γ, e, m)),
        )
    end
    return r
end

function _causality_residual(total::SelfEnergyFamily)
    return _maximum_block_residual(total.Σᴿ, _positive_imaginary_part)
end

function _spectral_sum_residual(problem::NEGFProblem, green::GreenState)
    Nb = problem.numerical.N_b
    I_b = Matrix{ComplexF64}(I, Nb, Nb)
    result = 0.0
    for m = 1:problem.numerical.N_k
        integral = zeros(ComplexF64, Nb, Nb)
        for e = 1:problem.numerical.N_E
            integral .+= problem.grids.wᴱ[e] / (2π) .* _matrix_block(green.A, e, m)
        end
        result = max(result, norm(integral - I_b) / sqrt(Nb))
    end
    return result
end

function _kernel_residual(K::AbstractArray{<:Number,6})
    Nk, Nk2, Nb, Nb2, Nb3, Nb4 = size(K)
    (Nk == Nk2 && Nb == Nb2 == Nb3 == Nb4) || return Inf, Inf
    rH = 0.0
    rPSD = 0.0
    compound = zeros(ComplexF64, Nb^2, Nb^2)
    for m = 1:Nk, mp = 1:Nk
        for a = 1:Nb, c = 1:Nb, b = 1:Nb, d = 1:Nb
            i = (a - 1) * Nb + c
            j = (b - 1) * Nb + d
            compound[i, j] = K[m, mp, a, c, d, b]
        end
        rH = max(rH, norm(compound - compound') / (norm(compound) + 1e-14))
        rPSD = max(rPSD, _negative_part(compound))
    end
    return rH, rPSD
end

function _nodal_shift_support_and_mass(W, ε, δ)
    support, mass = 0.0, 0.0
    N = length(ε)
    for i = 1:N
        target = ε[i] + δ
        inside = first(ε) <= target <= last(ε)
        left = inside ? clamp(searchsortedlast(ε, target), 1, N - 1) : 0
        for j = 1:N
            if !inside || (j != left && j != left + 1)
                support = max(support, abs(W[i, j]))
            end
        end
        mass = max(mass, abs(sum(view(W, i, :)) - (inside ? 1.0 : 0.0)))
    end
    return support, mass
end

"""Measure the weighted opposite-shift defect for the actual matrices.

This is a discretization diagnostic for every shift pair. Nodal interpolation
with zero exterior values is not quadrature adjoint at a finite-window edge;
only the conservative discretization promises roundoff-level equality.
"""
function _shift_quadrature_adjoint(plus, minus, weights)
    defect_squared, scale_squared = 0.0, 0.0
    for j in eachindex(weights), i in eachindex(weights)
        forward = weights[i] * plus[i, j]
        backward = weights[j] * minus[j, i]
        defect_squared += abs2(forward - backward)
        scale_squared += abs2(forward) + abs2(backward)
    end
    return iszero(scale_squared) ? 0.0 : sqrt(defect_squared / scale_squared)
end

"""Validate the support and integral identities of one declared shift pair.

Finite-volume support is determined by cell intersections. A shifted centre
outside the window does not imply an empty intersection. Both signs use the
same intersection geometry so their quadrature-adjoint relation is checked
without imposing the nodal interpolation boundary rule.
"""
function _shift_pair_validation(plus, minus, ε, weights, δ, discretization::Symbol)
    invalid = (
        wrap = Inf,
        nonnegative = Inf,
        row_mass = Inf,
        quadrature_adjoint = Inf,
        cell_weights = Inf,
    )
    N = length(ε)
    N >= 2 && size(plus) == size(minus) == (N, N) && length(weights) == N || return invalid
    all(isfinite, ε) &&
    all(>(0), diff(ε)) &&
    isfinite(δ) &&
    all(isfinite, weights) &&
    all(>(0), weights) &&
    all(isfinite, plus) &&
    all(isfinite, minus) || return invalid
    nonnegative = max(0.0, -minimum(plus), -minimum(minus))
    adjoint = _shift_quadrature_adjoint(plus, minus, weights)
    if discretization === :nodal_linear
        plus_support, plus_mass = _nodal_shift_support_and_mass(plus, ε, δ)
        minus_support, minus_mass = _nodal_shift_support_and_mass(minus, ε, -δ)
        return (
            wrap = max(plus_support, minus_support),
            nonnegative = nonnegative,
            row_mass = max(plus_mass, minus_mass),
            quadrature_adjoint = adjoint,
            cell_weights = 0.0,
        )
    elseif discretization !== :finite_volume_piecewise_constant
        return invalid
    end

    edges = vcat(first(ε), (ε[1:(end-1)] .+ ε[2:end]) ./ 2, last(ε))
    window_width = last(ε) - first(ε)
    support, mass = 0.0, 0.0
    for i = 1:N
        for j = 1:N
            overlap =
                max(0.0, min(edges[i+1] + δ, edges[j+1]) - max(edges[i] + δ, edges[j]))
            if iszero(overlap)
                support = max(support, abs(plus[i, j]), abs(minus[j, i]))
            end
        end
        expected_plus = max(0.0, min(edges[i+1] + δ, last(ε)) - max(edges[i] + δ, first(ε)))
        expected_minus =
            max(0.0, min(edges[i+1], last(ε) + δ) - max(edges[i], first(ε) + δ))
        # Integrate before normalizing by the energy window. Subtracting node
        # coordinates incurs coordinate roundoff, not one ulp of a tiny cell.
        mass = max(
            mass,
            abs(weights[i] * sum(view(plus, i, :)) - expected_plus) / window_width,
            abs(weights[i] * sum(view(minus, i, :)) - expected_minus) / window_width,
        )
    end
    return (
        wrap = support,
        nonnegative = nonnegative,
        row_mass = mass,
        quadrature_adjoint = adjoint,
        cell_weights = maximum(abs, weights .- diff(edges)) / window_width,
    )
end

function _antihermiticity_residual(X::AbstractArray{<:Number,4})
    return _maximum_block_residual(X, block -> norm(block + block') / (norm(block) + 1e-14))
end

function _selfenergy_identity_residual(total::SelfEnergyFamily)
    result = 0.0
    for e in axes(total.Σᴿ, 1), m in axes(total.Σᴿ, 2)
        ΣR = _matrix_block(total.Σᴿ, e, m)
        difference =
            ΣR - ΣR' - (_matrix_block(total.Σᵍ, e, m) - _matrix_block(total.Σˡ, e, m))
        denominator =
            norm(ΣR) +
            norm(_matrix_block(total.Σᵍ, e, m)) +
            norm(_matrix_block(total.Σˡ, e, m)) +
            1e-14
        result = max(result, norm(difference) / denominator)
    end
    return result
end

function _roundoff_metric(value::Complex, reversed::Complex)
    abs(value) > 1e-12 && return abs(value - reversed) / abs(value)
    return 100 * abs(value - reversed)
end

function _pointwise_current_imaginary(
    problem::NEGFProblem,
    family::SelfEnergyFamily,
    green::GreenState,
)
    maximum_imaginary = 0.0
    maximum_real = 0.0
    for e in axes(green.Gᴿ, 1), m in axes(green.Gᴿ, 2)
        value = _current_trace(family.Σᵍ, family.Σˡ, green.Gˡ, green.Gᵍ, e, m)
        maximum_imaginary = max(maximum_imaginary, abs(imag(value)))
        maximum_real = max(maximum_real, abs(real(value)))
    end
    return maximum_real > 1e-12 ? maximum_imaginary / maximum_real : 100 * maximum_imaginary
end

function _pointwise_collision_imaginary(family::SelfEnergyFamily, green::GreenState)
    maximum_imaginary = 0.0
    maximum_real = 0.0
    for e in axes(green.Gᴿ, 1), m in axes(green.Gᴿ, 2)
        incoming = tr(_matrix_block(family.Σˡ, e, m) * _matrix_block(green.Gᵍ, e, m))
        outgoing = tr(_matrix_block(family.Σᵍ, e, m) * _matrix_block(green.Gˡ, e, m))
        value = incoming - outgoing
        maximum_imaginary = max(maximum_imaginary, abs(imag(value)))
        maximum_real = max(maximum_real, abs(real(value)))
    end
    return maximum_real > 1e-12 ? maximum_imaginary / maximum_real : 100 * maximum_imaginary
end

function _edge_metrics(
    problem::NEGFProblem,
    green::GreenState,
    scattering::SelfEnergyFamily;
    energy_window_fraction::Real,
    momentum_window_fraction::Real,
)
    NE, Nk = size(green.Gᴿ, 1), size(green.Gᴿ, 2)
    Γsc = im .* (scattering.Σᵍ .- scattering.Σˡ)
    γall = 0.0
    γedge = 0.0
    Aall = 0.0
    Aedge = 0.0
    for e = 1:NE, m = 1:Nk
        γ = norm(_matrix_block(Γsc, e, m))
        a = abs(real(tr(_matrix_block(green.A, e, m))))
        γall = max(γall, γ)
        Aall = max(Aall, a)
        if e == 1 || e == NE
            γedge = max(γedge, γ)
            Aedge = max(Aedge, a)
        end
    end

    occupied = zeros(Float64, NE, Nk)
    for e = 1:NE, m = 1:Nk
        occupied[e, m] =
            problem.grids.wᴱ[e] *
            problem.grids.wᵏ[m] *
            real(-im * tr(_matrix_block(green.Gˡ, e, m)))
    end
    total = sum(occupied)
    if !isfinite(total) || total ≤ 0
        return (
            tail_low = Inf,
            tail_high = Inf,
            tail_k = Inf,
            edge_gamma = Inf,
            edge_spectral = Inf,
        )
    end
    nEedge = max(1, ceil(Int, energy_window_fraction * NE))
    nkedge = max(1, ceil(Int, momentum_window_fraction * Nk))
    return (
        tail_low = abs(sum(view(occupied, 1:nEedge, :))) / total,
        tail_high = abs(sum(view(occupied, (NE-nEedge+1):NE, :))) / total,
        tail_k = abs(sum(view(occupied, :, (Nk-nkedge+1):Nk))) / total,
        edge_gamma = γall == 0 ? 0.0 : γedge / γall,
        edge_spectral = Aall == 0 ? Inf : Aedge / Aall,
    )
end

const _STRUCTURAL_VALIDATION_LIMITS = Dict{Symbol,Float64}(
    :layer_sum => 1e-12,
    :donor_sum => 1e-12,
    :basis_orthogonality => 1e-12,
    :basis_eigen => 1e-12,
    :T_adjoint => 1e-12,
    :hamiltonian_hermiticity => 1e-12,
    :shift_wrap => 1e-14,
    :shift_nonnegative => 1e-14,
    :shift_row_mass => 1e-14,
    :shift_quadrature_adjoint => 1e-14,
    :shift_cell_weights => 1e-14,
    :shift_discretization => 0.0,
    :hilbert_trust_margin => 1e-14,
    :trusted_energy_nonempty => 0.0,
    :shape_contract => 0.0,
)
const _PZP_TRANSLATION_LIMIT = 0.1
const _PZP_TWO_PAIR_LIMIT = 5e-3
const _KERNEL_HERMITICITY_LIMIT = 1e-10
const _KERNEL_PSD_LIMIT = 1e-10
const _KERNEL_SCALE_LIMIT = 16eps(Float64)

"""
    validate_problem(problem)

Verify dimensions, donor sum, fixed-basis algebra, shift boundaries, projected
Hermiticity, and full-kernel covariance properties before an SCBA iteration is
allowed to start.  No clipping or silent repair is performed.

The shift operators carry their discretization in `NEGFProblem`. Both field
and LO pairs must be nonnegative and supported on their declared interpolation
nodes or intersecting control volumes. Conservative pairs must also satisfy
the quadrature-adjoint and cell-integral identities.

See [Verification](@ref theory-validation) and
[the array contract](@ref array-contracts).
"""
function validate_problem(problem::NEGFProblem)
    p, n, g = problem.physical, problem.numerical, problem.grids
    metrics = Dict{Symbol,Float64}()
    messages = String[]
    Lp = _metres(period_length(p))
    grid_length = sum(g.wˣ) * problem.scales.L₀_m
    metrics[:layer_sum] = abs(grid_length - Lp) / Lp
    target = _scaled_physics(p, problem.scales).Nᴰ²ᴰ
    metrics[:donor_sum] = abs(dot(g.wˣ, problem.profiles.Nᴰ) - target) / target
    metrics[:basis_orthogonality] = problem.basis.r_orth
    metrics[:basis_eigen] = problem.basis.r_eigen
    if problem.basis.localization in (:pzp,)
        metrics[:translation] = problem.basis.r_translation
        metrics[:translation_nearest] = problem.basis.r_translation_nearest
        metrics[:translation_two_pairs] = problem.basis.r_translation_two_pairs
    end
    metrics[:T_adjoint] =
        norm(problem.basis.T₋ - problem.basis.T₊') / (norm(problem.basis.T₊) + 1e-14)
    sp = _scaled_physics(p, problem.scales)
    h0 = project_hamiltonians(problem, zeros(n.N_z))
    metrics[:hamiltonian_hermiticity] = maximum(
        norm(Matrix(view(h0, m, :, :)) - Matrix(view(h0, m, :, :))') /
        (norm(Matrix(view(h0, m, :, :))) + 1e-14) for m = 1:n.N_k
    )
    discretization = problem.energy_shift_discretization
    metrics[:shift_discretization] =
        discretization in (:nodal_linear, :finite_volume_piecewise_constant) ? 0.0 : Inf
    field_shift =
        _shift_pair_validation(problem.W₊ᴱᵖ, problem.W₋ᴱᵖ, g.ε, g.wᴱ, sp.Eᵖ, discretization)
    phonon_shift = _shift_pair_validation(
        problem.W₊ᴸᴼ,
        problem.W₋ᴸᴼ,
        g.ε,
        g.wᴱ,
        sp.ħωᴸᴼ,
        discretization,
    )
    for name in keys(field_shift)
        metrics[Symbol(:shift_, name)] =
            max(getproperty(field_shift, name), getproperty(phonon_shift, name))
    end
    margin = scaled_value(n.M_E, problem.scales, :energy)
    required_margin = max(abs(sp.Eᵖ), sp.ħωᴸᴼ)
    metrics[:hilbert_trust_margin] =
        max(0.0, required_margin - margin) / max(required_margin, 1e-14)
    metrics[:trusted_energy_nonempty] = any(g.trusted_energy) ? 0.0 : Inf
    shapes_ok =
        size(problem.basis.Φ) == (n.N_z, n.N_b) &&
        all(
            size(W) == (n.N_E, n.N_E) for
            W in (problem.W₊ᴱᵖ, problem.W₋ᴱᵖ, problem.W₊ᴸᴼ, problem.W₋ᴸᴼ)
        ) &&
        length(g.trusted_energy) == n.N_E &&
        all(
            size(problem.kernels.K[name]) == (n.N_k, n.N_k, n.N_b, n.N_b, n.N_b, n.N_b) for
            name in problem.kernels.enabled
        )
    metrics[:shape_contract] = shapes_ok ? 0.0 : Inf
    for mechanism in problem.kernels.enabled
        rH, rPSD = _kernel_residual(problem.kernels.K[mechanism])
        metrics[Symbol(mechanism, :_hermiticity)] = rH
        metrics[Symbol(mechanism, :_PSD)] = rPSD
        metrics[Symbol(mechanism, :_kernel_scale)] =
            abs(maximum(abs, problem.kernels.K[mechanism]) - 1)
        q = problem.kernels.qᴷ[mechanism]
        metrics[Symbol(mechanism, :_kernel_scale_finite)] = isfinite(q) && q > 0 ? 0.0 : Inf
    end
    limits = copy(_STRUCTURAL_VALIDATION_LIMITS)
    if discretization === :nodal_linear
        # Keep the measured defect in the report. The exact weighted-adjoint
        # structural contract belongs to finite volumes, not nodal sampling.
        # Scientific collision and power gates are unchanged and still apply.
        delete!(limits, :shift_quadrature_adjoint)
        metrics[:shift_quadrature_adjoint] >
        _STRUCTURAL_VALIDATION_LIMITS[:shift_quadrature_adjoint] && push!(
            messages,
            "nodal shift quadrature-adjoint defect=$(metrics[:shift_quadrature_adjoint]); measured discretization error, assess collision and grid convergence",
        )
    end
    if problem.basis.localization in (:pzp,)
        limits[:translation] = _PZP_TRANSLATION_LIMIT
        limits[:translation_nearest] = _PZP_TRANSLATION_LIMIT
        limits[:translation_two_pairs] =
            n.P_basis ≥ 3 ? _PZP_TWO_PAIR_LIMIT : _PZP_TRANSLATION_LIMIT
    end
    passed = true
    for (name, limit) in limits
        if !isfinite(metrics[name]) || metrics[name] > limit
            passed = false
            push!(messages, "$name=$(metrics[name]) exceeds $limit")
        end
    end
    for mechanism in problem.kernels.enabled
        for (suffix, limit) in (
            (:_hermiticity, _KERNEL_HERMITICITY_LIMIT),
            (:_PSD, _KERNEL_PSD_LIMIT),
            (:_kernel_scale, _KERNEL_SCALE_LIMIT),
            (:_kernel_scale_finite, 0.0),
        )
            name = Symbol(mechanism, suffix)
            if !isfinite(metrics[name]) || metrics[name] > limit
                passed = false
                push!(messages, "$name=$(metrics[name]) exceeds $limit")
            end
        end
    end
    for (name, value) in metrics
        if !isfinite(value) && !any(occursin(String(name), message) for message in messages)
            passed = false
            push!(messages, "$name is non-finite")
        end
    end
    return ConvergenceReport(passed, metrics, messages)
end

"""
    validate(problem_or_solution)

Run the non-mutating static or converged-state validation suite and return a
[`ConvergenceReport`](@ref).  No matrix is symmetrized, clipped, or repaired.

See [Verification](@ref theory-validation).
"""
validate(problem::NEGFProblem) = validate_problem(problem)

"""
Run algebraic, conservation, current, power, and spectral checks.

See [Verification](@ref theory-validation), including the complete list of
acceptance residuals and their physical interpretation.
"""
function validate_solution(
    solution::NEGFSolution;
    candidate::Union{Nothing,FixedPointAuditCandidate} = nothing,
    audit_failure::Union{Nothing,AbstractString} = nothing,
)
    problem = solution.problem
    scba = solution.scba
    green = scba.green
    total_sc = _sum_selfenergies(scba.scattering, size(green.Gᴿ))
    total = SelfEnergyFamily(
        total_sc.Σᴿ + scba.embedding.Σᴿ,
        total_sc.Σˡ + scba.embedding.Σˡ,
        total_sc.Σᵍ + scba.embedding.Σᵍ,
    )
    metrics = Dict{Symbol,Float64}()
    if isempty(scba.history)
        for name in
            (:r_D, :r_A, :r_K, :r_Σ, :r_λ, :r_roundoff, :r_Jchange_scba, :r_population_scba)
            metrics[name] = Inf
        end
    else
        last_inner = last(scba.history)
        metrics[:r_D] = last_inner.r_D
        metrics[:r_A] = last_inner.r_A
        metrics[:r_K] = last_inner.r_K
        metrics[:r_Σ] = last_inner.r_Σ
        metrics[:r_λ] = last_inner.r_λ
        metrics[:r_roundoff] = last_inner.r_roundoff
        metrics[:r_Jchange_scba] = last_inner.r_Jchange
        metrics[:r_population_scba] = last_inner.r_population
    end
    # Never certify from stored iteration metrics alone. Production callers
    # supply their explicitly selected map; direct reference calls rebuild it.
    if audit_failure === nothing && scba.quality !== :invalid
        candidate === nothing && (candidate = _audit_candidate_for_solution(solution))
        merge!(metrics, fresh_fixed_point_metrics(problem, scba, solution.Uᴴ, candidate))
    else
        # Preserve the aligned state when the selected nonlinear model is
        # undefined on it (e.g. a negative hot-phonon event rate). Missing fresh
        # fixed-point checks are failures, never stale green certificates.
        metrics[:r_D] = _dyson_residual(
            problem.grids.ε,
            project_hamiltonians(problem, solution.Uᴴ),
            total.Σᴿ,
            green.Gᴿ,
        )
        metrics[:r_A] = _spectral_residual(green, total)
        metrics[:r_K] = Inf
        metrics[:r_Σ] = Inf
        metrics[:r_λ] = Inf
    end
    metrics[:r_PSD] = _green_positivity(green, total)
    metrics[:r_caus] = _causality_residual(total)
    metrics[:r_caus_sc] = _causality_residual(total_sc)
    metrics[:r_caus_embedding] = _causality_residual(scba.embedding)
    metrics[:r_selfenergy_identity] = _selfenergy_identity_residual(total)
    metrics[:r_selfenergy_identity_sc] = _selfenergy_identity_residual(total_sc)
    metrics[:r_selfenergy_identity_embedding] =
        _selfenergy_identity_residual(scba.embedding)
    metrics[:r_antihermitian_G] =
        max(_antihermiticity_residual(green.Gˡ), _antihermiticity_residual(green.Gᵍ))
    metrics[:r_antihermitian_Σ] =
        max(_antihermiticity_residual(total.Σˡ), _antihermiticity_residual(total.Σᵍ))
    metrics[:r_antihermitian_Σsc] =
        max(_antihermiticity_residual(total_sc.Σˡ), _antihermiticity_residual(total_sc.Σᵍ))
    metrics[:r_antihermitian_Σembedding] = max(
        _antihermiticity_residual(scba.embedding.Σˡ),
        _antihermiticity_residual(scba.embedding.Σᵍ),
    )
    metrics[:r_sum] = _spectral_sum_residual(problem, green)
    for (name, family) in scba.scattering
        balance = collision_balance(problem, family, green)
        metrics[Symbol(:r_C_, name)] = balance.residual
        metrics[Symbol(:r_ImC_, name)] = _pointwise_collision_imaginary(family, green)
        metrics[Symbol(:r_roundC_, name)] = balance.roundoff
    end
    if :electron_electron in problem.kernels.enabled
        ee = electron_electron_collision_diagnostics(
            problem,
            green,
            scba.scattering[:electron_electron],
        )
        metrics[:r_energy_electron_electron] = ee.energy_residual
        metrics[:electron_electron_energy_transfer_W_m2] = ee.electron_energy_rate_W_m2
        metrics[:electron_electron_particle_rate_m2_s] = ee.electron_number_rate_m2_s
        # A prescribed single-q plasmon bath can exchange electronic energy.
        # Its signed power remains part of the total JV + Pscattering balance;
        # vanishing electronic power alone is not its conservation law.
        metrics[:ee_energy_conservation_applicable] =
            problem.models.electron_electron.mode === :sppa_single_q ? 0.0 : 1.0
    end
    if problem.scattering.LO && problem.models.lo_population === :rate_balance
        metrics[:r_hot_lo_power] = try
            _lo_balance_diagnostics(problem, scba).collision_balance_residual
        catch error
            error isa DomainError || rethrow()
            Inf
        end
    end
    Φplus = _boundary_flux_complex(problem, scba.embedding_plus, green)
    Φminus = _boundary_flux_complex(problem, scba.embedding_minus, green)
    metrics[:r_J] =
        abs(real(Φplus + Φminus)) / max(abs(real(Φplus)), abs(real(Φminus)), 1e-12)
    metrics[:r_ImJ] = max(
        _pointwise_current_imaginary(problem, scba.embedding_plus, green),
        _pointwise_current_imaginary(problem, scba.embedding_minus, green),
    )
    Φplus_reverse =
        _boundary_flux_complex(problem, scba.embedding_plus, green; reverse_order = true)
    Φminus_reverse =
        _boundary_flux_complex(problem, scba.embedding_minus, green; reverse_order = true)
    metrics[:r_roundJ] = max(
        _roundoff_metric(Φplus, Φplus_reverse),
        _roundoff_metric(Φminus, Φminus_reverse),
    )
    pbalance = power_balance(problem, scba.scattering, scba.embedding_plus, green)
    metrics[:r_power] = pbalance.residual
    metrics[:r_ImP] = 100 * pbalance.imaginary
    metrics[:r_roundP] = pbalance.roundoff

    Up, ζ, rP, rneutral = solve_periodic_poisson(problem, solution.n)
    metrics[:r_P] = rP
    metrics[:r_U] = maximum(abs.(Up .- solution.Uᴴ))
    metrics[:r_neutral] = rneutral
    metrics[:r_ζ] = abs(ζ)
    if isempty(solution.outer_history)
        metrics[:r_n] = Inf
        metrics[:r_Jchange] = Inf
        metrics[:r_population] = Inf
    else
        last_outer = last(solution.outer_history)
        density_normalizer =
            _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ /
            sqrt(_scaled_physics(problem.physical, problem.scales).Lp)
        final_density_change =
            sqrt(sum(problem.grids.wˣ .* (solution.n .- last_outer.density) .^ 2)) /
            density_normalizer
        metrics[:r_n] = max(last_outer.r_n, final_density_change)
        J₀ = Float64(ustrip(u"A/m^2", problem.scales.J₀))
        Jfinal = real(Φplus) * J₀
        final_Jchange =
            abs(Jfinal - last_outer.J) / max(abs(Jfinal), abs(last_outer.J), 1e-12)
        target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
        final_populations =
            real.(diag(_sheet_density_matrix_bar(problem, green.Gˡ))) ./ target
        final_population_change =
            norm(final_populations - last_outer.populations) /
            (norm(final_populations) + 1e-14)
        metrics[:r_Jchange] = max(last_outer.r_Jchange, final_Jchange)
        metrics[:r_population] = max(last_outer.r_population, final_population_change)
    end
    edges = _edge_metrics(
        problem,
        green,
        total_sc;
        energy_window_fraction = solution.options.energy_tail_window_fraction,
        momentum_window_fraction = solution.options.momentum_tail_window_fraction,
    )
    metrics[:r_tail_low] = edges.tail_low
    metrics[:r_tail_high] = edges.tail_high
    metrics[:r_tail_k] = edges.tail_k
    metrics[:r_edge_gamma] = edges.edge_gamma
    metrics[:r_edge_spectral] = edges.edge_spectral
    report = _validate_solution_metrics(solution, metrics)
    if audit_failure !== nothing || scba.quality === :invalid
        push!(
            report.messages,
            "Fresh physical map audit unavailable: " *
            (audit_failure === nothing ? String(scba.status) : String(audit_failure)),
        )
    end
    return report
end

"""Declared final-state thresholds, shared by validation and light diagnostics."""
function _solution_validation_limits(t::SolverTolerances, metrics)
    limits = Dict(
        :r_D => t.r_D,
        :r_A => t.r_A,
        :r_K => t.r_K,
        :r_Σ => t.r_Σ,
        :r_λ => t.r_λ,
        :r_roundoff => t.r_roundoff,
        :r_Jchange_scba => t.r_obs,
        :r_population_scba => t.r_obs,
        :r_PSD => t.r_PSD,
        :r_caus => t.r_caus,
        :r_caus_sc => t.r_caus,
        :r_caus_embedding => t.r_caus,
        :r_selfenergy_identity => t.r_caus,
        :r_selfenergy_identity_sc => t.r_caus,
        :r_selfenergy_identity_embedding => t.r_caus,
        :r_antihermitian_G => t.r_caus,
        :r_antihermitian_Σ => t.r_caus,
        :r_antihermitian_Σsc => t.r_caus,
        :r_antihermitian_Σembedding => t.r_caus,
        :r_sum => t.r_sum,
        :r_J => t.r_J,
        :r_ImJ => t.r_imag,
        :r_roundJ => t.r_roundoff,
        :r_power => t.r_power,
        :r_ImP => t.r_imag,
        :r_roundP => t.r_roundoff,
        :r_P => t.r_P,
        :r_U => t.r_U,
        :r_n => t.r_n,
        :r_neutral => t.r_neutral,
        :r_ζ => t.r_ζ,
        :r_Jchange => t.r_obs,
        :r_population => t.r_obs,
        :r_tail_low => t.r_tail,
        :r_tail_high => t.r_tail,
        :r_tail_k => t.r_tail,
        :r_edge_gamma => t.r_edge,
        :r_edge_spectral => t.r_edge,
    )
    for name in keys(metrics)
        label = String(name)
        startswith(label, "r_C_") && (limits[name] = t.r_C)
        startswith(label, "r_ImC_") && (limits[name] = t.r_imag)
        startswith(label, "r_roundC_") && (limits[name] = t.r_roundoff)
    end
    haskey(metrics, :r_energy_electron_electron) &&
        get(metrics, :ee_energy_conservation_applicable, 1.0) == 1.0 &&
        (limits[:r_energy_electron_electron] = t.r_power)
    haskey(metrics, :r_hot_lo_power) && (limits[:r_hot_lo_power] = t.r_power)
    haskey(metrics, :r_lambda_candidate) && (limits[:r_lambda_candidate] = t.r_λ)
    # Candidate equations belong to the nonlinear map, distinct from the
    # machine-precision linear solve with the stored self-energy.
    haskey(metrics, :r_D_candidate) && (limits[:r_D_candidate] = t.r_K)
    haskey(metrics, :r_A_candidate) && (limits[:r_A_candidate] = t.r_K)
    return limits
end

"""Reassess the same final state; approximate acceptance never changes its SCBA label."""
function _validate_solution_metrics(
    solution::NEGFSolution,
    metrics;
    accept_approximate_scba::Bool = false,
)
    limits = _solution_validation_limits(solution.options.tolerances, metrics)
    accepted =
        solution.scba.converged || (
            accept_approximate_scba &&
            solution.options.convergence.diagnostic_quality.enabled &&
            scba_accepted(solution.scba)
        )
    messages =
        accepted ? String[] :
        ["SCBA state is not accepted for this quality band: $(solution.scba.status)"]
    passed = accepted
    for name in sort!(collect(keys(limits)); by = String)
        limit = limits[name]
        value = get(metrics, name, Inf)
        if !isfinite(value) || value > limit
            passed = false
            push!(messages, "$name=$value exceeds $limit")
        end
    end
    return ConvergenceReport(passed, metrics, messages)
end

validate(solution::NEGFSolution) = validate_solution(solution)
