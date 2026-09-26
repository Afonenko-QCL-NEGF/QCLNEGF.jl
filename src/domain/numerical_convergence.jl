"""Validate the solver options shared by reference and optimized numerics."""
function _check_options(options::SolverOptions)
    0 < options.α_Σ ≤ 1 || throw(ArgumentError("α_Σ must lie in (0,1]"))
    0 < options.α_P ≤ 1 || throw(ArgumentError("α_P must lie in (0,1]"))
    options.max_scba ≥ 1 || throw(ArgumentError("max_scba must be positive"))
    options.max_poisson ≥ 1 || throw(ArgumentError("max_poisson must be positive"))
    0 < options.energy_tail_window_fraction < 1 ||
        throw(ArgumentError("energy_tail_window_fraction must lie in (0,1)"))
    0 < options.momentum_tail_window_fraction < 1 ||
        throw(ArgumentError("momentum_tail_window_fraction must lie in (0,1)"))
    for name in fieldnames(SolverTolerances)
        value = getfield(options.tolerances, name)
        isfinite(value) && value ≥ 0 ||
            throw(ArgumentError("tolerance $name must be finite and nonnegative"))
    end
    policy = options.convergence
    policy.mode in (:strict_fail_fast, :research_continue, :adaptive_working) ||
        throw(ArgumentError("unknown continuation policy $(policy.mode)"))
    policy.minimum_scba_iterations >= 1 ||
        throw(ArgumentError("minimum_scba_iterations must be positive"))
    policy.required_consecutive_scba_passes >= 1 ||
        throw(ArgumentError("required_consecutive_scba_passes must be positive"))
    policy.minimum_poisson_iterations >= 1 ||
        throw(ArgumentError("minimum_poisson_iterations must be positive"))
    policy.required_consecutive_poisson_passes >= 1 ||
        throw(ArgumentError("required_consecutive_poisson_passes must be positive"))
    policy.stagnation_window == 0 ||
        policy.stagnation_window >= 4 ||
        throw(ArgumentError("stagnation_window must be zero or at least four"))
    0 <= policy.stagnation_relative_improvement < 1 ||
        throw(ArgumentError("stagnation_relative_improvement must lie in [0,1)"))
    policy.stagnation_window == 0 &&
        policy.stagnation_relative_improvement != 0 &&
        throw(
            ArgumentError(
                "stagnation_relative_improvement requires a nonzero stagnation_window",
            ),
        )
    quality = policy.diagnostic_quality
    for (name, value) in (
        (:keldysh_threshold, quality.keldysh_threshold),
        (:self_energy_threshold, quality.self_energy_threshold),
        (:normalization_threshold, quality.normalization_threshold),
        (:target_fixed_point_threshold, quality.target_fixed_point_threshold),
        (:observable_threshold, quality.observable_threshold),
        (:positivity_threshold, quality.positivity_threshold),
        (:causality_threshold, quality.causality_threshold),
        (:outer_threshold, quality.outer_threshold),
    )
        isfinite(value) && value > 0 ||
            throw(ArgumentError("$name must be finite and positive"))
    end
    quality.required_consecutive_passes >= 3 || throw(
        ArgumentError(
            "diagnostic quality required_consecutive_passes must be at least three",
        ),
    )
    quality.coarse_wait_iterations >= 0 ||
        throw(ArgumentError("coarse_wait_iterations cannot be negative"))
    quality.strict_attempt_iterations >= 0 ||
        throw(ArgumentError("strict_attempt_iterations cannot be negative"))
    quality.target_fixed_point_threshold <=
    min(quality.keldysh_threshold, quality.self_energy_threshold) || throw(
        ArgumentError(
            "target_fixed_point_threshold cannot exceed the coarse fixed-point limits",
        ),
    )
    return options
end

const _SCBA_RESULT_QUALITIES =
    (:strictly_converged, :approximate_fixed_point, :unresolved, :invalid)

"""Validate the shared outer-status and SCBA-quality classification."""
function _check_scba_result_classification(converged::Bool, status::Symbol, quality::Symbol)
    quality in _SCBA_RESULT_QUALITIES ||
        throw(ArgumentError("unsupported SCBA quality $quality"))
    converged == (status === :converged) ||
        throw(ArgumentError("converged must be true if and only if status is :converged"))
    !converged ||
        quality === :strictly_converged ||
        throw(ArgumentError("a converged result requires strictly_converged SCBA quality"))
    return quality
end

"""Validate the stricter classification of one fixed-Hartree SCBA state."""
function _check_inner_scba_result_classification(
    converged::Bool,
    status::Symbol,
    quality::Symbol,
)
    _check_scba_result_classification(converged, status, quality)
    converged ||
        quality !== :strictly_converged ||
        throw(ArgumentError("a non-converged SCBA state cannot be strictly_converged"))
    return quality
end

"""Result of evaluating every required gate at one solver iteration."""
struct ConvergenceAssessment
    eligible::Bool
    passed::Bool
    limiting_metric::Symbol
    limiting_ratio::Float64
    failed_metrics::Vector{Symbol}
end

const _SCBA_CONVERGENCE_METRICS = (
    (:r_D, :r_D),
    (:r_A, :r_A),
    (:r_K, :r_K),
    (:r_Σ, :r_Σ),
    (:r_λ, :r_λ),
    (:r_PSD, :r_PSD),
    (:r_caus, :r_caus),
    (:r_roundoff, :r_roundoff),
    (:r_Jchange, :r_obs),
    (:r_population, :r_obs),
)

const _POISSON_CONVERGENCE_METRICS = (
    (:r_P, :r_P),
    (:r_U, :r_U),
    (:r_n, :r_n),
    (:r_neutral, :r_neutral),
    (:r_J, :r_J),
    (:r_Jchange, :r_obs),
    (:r_population, :r_obs),
    (:ζ, :r_ζ),
)

@inline function _gate_value(iteration, metric::Symbol)
    value = Float64(getproperty(iteration, metric))
    return metric === :ζ ? abs(value) : value
end

@inline function _gate_ratio(value::Float64, threshold::Float64)
    isfinite(value) || return Inf
    threshold > 0 || return value < threshold ? 0.0 : Inf
    return value / threshold
end

function _assess_convergence(
    iteration,
    tolerances::SolverTolerances,
    metrics,
    minimum_iteration::Integer,
)
    limiting_metric = first(metrics)[1]
    limiting_ratio = -Inf
    failed = Symbol[]
    for (metric, tolerance_field) in metrics
        value = _gate_value(iteration, metric)
        threshold = getfield(tolerances, tolerance_field)
        ratio = _gate_ratio(value, threshold)
        if ratio > limiting_ratio
            limiting_metric = metric
            limiting_ratio = ratio
        end
        value < threshold || push!(failed, metric)
    end
    eligible =
        getproperty(iteration, iteration isa SCBAIteration ? :ν : :μ) >= minimum_iteration
    return ConvergenceAssessment(
        eligible,
        eligible && isempty(failed),
        limiting_metric,
        limiting_ratio,
        failed,
    )
end

"""Evaluate the complete inner fixed-point acceptance contract."""
scba_convergence_assessment(
    iteration::SCBAIteration,
    tolerances::SolverTolerances,
    policy::ConvergencePolicy,
) = _assess_convergence(
    iteration,
    tolerances,
    _SCBA_CONVERGENCE_METRICS,
    policy.minimum_scba_iterations,
)

"""Evaluate the complete outer Poisson fixed-point acceptance contract."""
poisson_convergence_assessment(
    iteration::OuterIteration,
    tolerances::SolverTolerances,
    policy::ConvergencePolicy,
) = _assess_convergence(
    iteration,
    tolerances,
    _POISSON_CONVERGENCE_METRICS,
    policy.minimum_poisson_iterations,
)

"""
Evaluate the optional exploratory SCBA quality band.

Fixed-point, normalization and declared physical stability gates use their
separate approximate limits. Dyson, spectral identity and roundoff retain
the strict algebraic tolerances. All measurements belong to the current state.
"""
function scba_diagnostic_quality_assessment(
    iteration::SCBAIteration,
    tolerances::SolverTolerances,
    policy::ConvergencePolicy,
)
    quality = policy.diagnostic_quality
    quality.enabled || return ConvergenceAssessment(
        false,
        false,
        :diagnostic_quality_disabled,
        Inf,
        Symbol[],
    )
    limiting_metric = first(_SCBA_CONVERGENCE_METRICS)[1]
    limiting_ratio = -Inf
    failed = Symbol[]
    for (metric, tolerance_field) in _SCBA_CONVERGENCE_METRICS
        threshold =
            metric === :r_K ? quality.keldysh_threshold :
            metric === :r_Σ ? quality.self_energy_threshold :
            metric === :r_λ ? quality.normalization_threshold :
            metric === :r_PSD ? quality.positivity_threshold :
            metric === :r_caus ? quality.causality_threshold :
            metric in (:r_Jchange, :r_population) ? quality.observable_threshold :
            getfield(tolerances, tolerance_field)
        value = _gate_value(iteration, metric)
        ratio = _gate_ratio(value, threshold)
        if ratio > limiting_ratio
            limiting_metric = metric
            limiting_ratio = ratio
        end
        value < threshold || push!(failed, metric)
    end
    eligible = iteration.ν >= policy.minimum_scba_iterations
    return ConvergenceAssessment(
        eligible,
        eligible && isempty(failed),
        limiting_metric,
        limiting_ratio,
        failed,
    )
end

"""Advance the approximate-band diagnostic streak and derive its quality label."""
function _scba_diagnostic_quality_state(
    previous_streak::Integer,
    assessment::ConvergenceAssessment,
    policy::ConvergencePolicy,
)
    streak =
        policy.diagnostic_quality.enabled && assessment.passed ? previous_streak + 1 : 0
    quality =
        streak >= policy.diagnostic_quality.required_consecutive_passes ?
        :approximate_fixed_point : :unresolved
    return (; streak, quality)
end

# Physical admissibility is independent of progress of the nonlinear map.
# In particular, the normalized PSD defect is bounded near one and can stay
# saturated even while the raw fixed-point residual decreases substantially.
const _SCBA_FIXED_POINT_METRICS = (
    (:r_K, :r_K),
    (:r_Σ, :r_Σ),
    (:r_λ, :r_λ),
    (:r_Jchange, :r_obs),
    (:r_population, :r_obs),
)
const _SCBA_ADMISSIBILITY_METRICS = (
    (:r_D, :r_D),
    (:r_A, :r_A),
    (:r_PSD, :r_PSD),
    (:r_caus, :r_caus),
    (:r_roundoff, :r_roundoff),
)

_scba_fixed_point_assessment(row, tolerances, policy) = _assess_convergence(
    row,
    tolerances,
    _SCBA_FIXED_POINT_METRICS,
    policy.minimum_scba_iterations,
)

function _scba_window_progress(history, window, score)
    first_index = length(history) - window + 1
    split_index = first_index + window ÷ 2 - 1
    old_best, new_best = Inf, Inf
    for index = first_index:length(history)
        value = score(history[index])
        # A nonfinite window cannot certify either stagnation or validity.
        isfinite(value) || return nothing
        if index <= split_index
            old_best = min(old_best, value)
        else
            new_best = min(new_best, value)
        end
    end
    improvement = max(0.0, (old_best-new_best) / max(old_best, eps(Float64)))
    return (; old_best, new_best, improvement)
end

"""Diagnose bounded failure without using physical quality as a progress score.

`:stagnated` requires an unaccepted fixed-point score to stop improving.
`:quality_blocked` requires fixed-point gates to pass throughout the window
while every remaining inadmissibility gate fails to improve. Neither status
accepts the state or relaxes any strict or approximate tolerance.
"""
function _scba_stopping_reason(
    history::AbstractVector{SCBAIteration},
    tolerances::SolverTolerances,
    policy::ConvergencePolicy,
)
    policy.mode === :research_continue && return nothing
    window = policy.stagnation_window
    window > 0 && length(history) >= window || return nothing
    progress = _scba_window_progress(
        history,
        window,
        row -> _scba_fixed_point_assessment(row, tolerances, policy).limiting_ratio,
    )
    progress === nothing && return nothing
    if progress.new_best >= 1
        return progress.improvement < policy.stagnation_relative_improvement ? :stagnated :
               nothing
    end
    # A temporarily small score is not a converged fixed-point sequence.
    all(
        row -> _scba_fixed_point_assessment(row, tolerances, policy).passed,
        @view history[(end-window+1):end]
    ) || return nothing
    remaining = [
        (metric, threshold) for (metric, threshold) in _SCBA_ADMISSIBILITY_METRICS if
        !(_gate_value(last(history), metric) < getfield(tolerances, threshold))
    ]
    isempty(remaining) && return nothing
    for (metric, threshold) in remaining
        quality_progress = _scba_window_progress(
            history,
            window,
            row ->
                _gate_ratio(_gate_value(row, metric), getfield(tolerances, threshold)),
        )
        quality_progress === nothing && return nothing
        quality_progress.new_best >= 1 || return nothing
        quality_progress.improvement < policy.stagnation_relative_improvement ||
            return nothing
    end
    return :quality_blocked
end

"""Whether the SCBA fixed-point iteration, independently of physical gates, stagnates."""
scba_convergence_stagnated(
    history::AbstractVector{SCBAIteration},
    tolerances::SolverTolerances,
    policy::ConvergencePolicy,
) = _scba_stopping_reason(history, tolerances, policy) === :stagnated


function _trailing_pass_count(history, passed::Function)
    count = 0
    for row in Iterators.reverse(history)
        passed(row) || break
        count += 1
    end
    return count
end

"""Bounded approximate acceptance, based on fresh raw residuals before mixing.

The complete last confirmation window must pass. The target is attempted
first; a continuously valid coarse band is accepted after a bounded wait.
No change of mixing strength can change these raw residual gates.
"""
function scba_approximate_acceptance(
    history::AbstractVector{SCBAIteration},
    options::SolverOptions,
)
    policy = options.convergence
    q = policy.diagnostic_quality
    policy.mode === :adaptive_working || return false
    q.enabled || return false
    n = length(history)
    n >= q.required_consecutive_passes || return false
    coarse_streak = 0
    target_streak = 0
    for row in Iterators.reverse(history)
        a = scba_diagnostic_quality_assessment(row, options.tolerances, policy)
        a.passed || break
        coarse_streak += 1
        if target_streak == coarse_streak - 1 &&
           max(row.r_K, row.r_Σ) < q.target_fixed_point_threshold
            target_streak += 1
        end
    end
    target_ready =
        target_streak >= q.required_consecutive_passes + q.strict_attempt_iterations
    coarse_ready = coarse_streak >= q.required_consecutive_passes + q.coarse_wait_iterations
    # A budget edge can use a fully confirmed band, never a single crossing.
    exhausted = !isempty(history) && last(history).ν >= options.max_scba
    stabilized = scba_convergence_stagnated(history, options.tolerances, policy)
    return target_ready ||
           coarse_ready ||
           ((exhausted || stabilized) && coarse_streak >= q.required_consecutive_passes)
end

"""Whether an inner result permits evaluating the outer fixed-point map."""
scba_accepted(scba::SCBAResult) =
    scba.converged ||
    (scba.status === :research_continue && scba.quality === :unresolved) ||
    (scba.status === :approximate && scba.quality === :approximate_fixed_point)

"""Stable public quality vocabulary, independent of process completion."""
solution_quality(solution::NEGFSolution) =
    solution.converged ? :strict :
    solution.scba.quality === :invalid ? :invalid :
    solution.status === :approximate ? :approximate : :unconverged

function _approximate_options(options::SolverOptions)
    q = options.convergence.diagnostic_quality
    t = options.tolerances
    changes = (
        r_K = q.keldysh_threshold,
        r_Σ = q.self_energy_threshold,
        r_λ = q.normalization_threshold,
        r_PSD = q.positivity_threshold,
        r_caus = q.causality_threshold,
        r_obs = q.observable_threshold,
        r_U = max(t.r_U, q.outer_threshold),
        r_n = max(t.r_n, q.outer_threshold),
    )
    tolerances = SolverTolerances(;
        (
            name => (
                hasproperty(changes, name) ? getproperty(changes, name) : getfield(t, name)
            ) for name in fieldnames(SolverTolerances)
        )...,
    )
    return SolverOptions(;
        (
            name => (name === :tolerances ? tolerances : getfield(options, name)) for
            name in fieldnames(SolverOptions)
        )...,
    )
end

"""Choose bounded inexact inner tolerances from completed raw outer residuals.

Only `adaptive_working` with an explicitly enabled diagnostic band uses this
schedule. The outer solve always retains its strict tolerances. The best
finite raw Hartree/density residual tightens the inner band monotonically;
within ten strict outer tolerances the approximate return is disabled.
This does not change the physical gates or the final certification thresholds.
"""
function inner_working_policy(
    options::SolverOptions,
    outer_history::AbstractVector{OuterIteration},
)
    policy = options.convergence
    q = policy.diagnostic_quality
    if policy.mode !== :adaptive_working || !q.enabled
        return (;
            options,
            stage = :strict,
            forcing_threshold = nothing,
            outer_residual = nothing,
        )
    end
    finite_rows = filter(row -> isfinite(row.r_U) && isfinite(row.r_n), outer_history)
    best =
        isempty(finite_rows) ? Inf : minimum(max(row.r_U, row.r_n) for row in finite_rows)
    strict_ready = any(
        row ->
            row.r_U <= 10options.tolerances.r_U && row.r_n <= 10options.tolerances.r_n,
        finite_rows,
    )
    forcing = isfinite(best) ? 0.1best : max(q.keldysh_threshold, q.self_energy_threshold)
    t = options.tolerances
    changes = if strict_ready
        (; enabled = false)
    else
        keldysh = max(t.r_K, min(q.keldysh_threshold, forcing))
        self_energy = max(t.r_Σ, min(q.self_energy_threshold, forcing))
        (;
            keldysh_threshold = keldysh,
            self_energy_threshold = self_energy,
            normalization_threshold = max(t.r_λ, min(q.normalization_threshold, forcing)),
            observable_threshold = max(t.r_obs, min(q.observable_threshold, forcing)),
            target_fixed_point_threshold = min(
                q.target_fixed_point_threshold,
                keldysh,
                self_energy,
            ),
        )
    end
    quality = DiagnosticQualityPolicy(;
        (
            name =>
                hasproperty(changes, name) ? getproperty(changes, name) : getfield(q, name) for name in fieldnames(DiagnosticQualityPolicy)
        )...,
    )
    effective_policy = ConvergencePolicy(;
        (
            name => name === :diagnostic_quality ? quality : getfield(policy, name) for
            name in fieldnames(ConvergencePolicy)
        )...,
    )
    effective_options = SolverOptions(;
        (
            name => name === :convergence ? effective_policy : getfield(options, name)
            for name in fieldnames(SolverOptions)
        )...,
    )
    return (;
        options = effective_options,
        stage = strict_ready ? :strict_final : :inexact_inner,
        forcing_threshold = strict_ready ? nothing : forcing,
        outer_residual = isfinite(best) ? best : nothing,
    )
end

function _approximate_warning(
    scba::SCBAResult,
    options::SolverOptions,
    outer_iteration::Integer,
)
    q = options.convergence.diagnostic_quality
    row = last(scba.history)
    return Dict{String,Any}(
        "code"=>"SCBA_APPROXIMATE_ACCEPTED",
        "scope"=>"inner_scba",
        "outer_iteration"=>Int(outer_iteration),
        "iteration"=>row.ν,
        "message"=>"Approximate fixed-Hartree state accepted; strict convergence is not certified.",
        "thresholds"=>Dict(
            "target"=>q.target_fixed_point_threshold,
            "keldysh"=>q.keldysh_threshold,
            "self_energy"=>q.self_energy_threshold,
            "normalization"=>q.normalization_threshold,
            "positivity"=>q.positivity_threshold,
            "causality"=>q.causality_threshold,
            "observables"=>q.observable_threshold,
            "confirmation_passes"=>q.required_consecutive_passes,
        ),
        "metrics"=>Dict(
            "r_K"=>row.r_K,
            "r_Sigma"=>row.r_Σ,
            "lambda_minus_one"=>row.λ-1,
            "r_PSD"=>row.r_PSD,
            "r_caus"=>row.r_caus,
            "r_Jchange"=>row.r_Jchange,
            "r_population"=>row.r_population,
        ),
    )
end
