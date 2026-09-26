"""One reproducible summary of a scalar residual, without physical gate mixing."""
struct ResidualTrend
    initial::Float64
    final::Float64
    initial_window::Float64
    final_window::Float64
    best_window::Float64
    best_value::Float64
    reduction::Float64
    regression::Float64
    oscillation::Float64
    last_change::Float64
    promising::Bool
end

"""Permission for an outer step is independent of strict scientific acceptance."""
struct ResearchContinuationDecision
    action::Symbol
    reason::Symbol
    iteration::Int
    budget::Int
    usable::Bool
    selected_state::Symbol
    keldysh::ResidualTrend
    self_energy::ResidualTrend
    normalization::ResidualTrend
end

function residual_trend(values::AbstractVector{<:Real}, policy::ResearchTrendPolicy)
    invalid = ResidualTrend(NaN, NaN, NaN, NaN, NaN, NaN, 0.0, Inf, Inf, Inf, false)
    isempty(values) && return invalid
    all(x -> isfinite(x) && x >= 0, values) || return invalid
    n = length(values)
    width = min(policy.window, n)
    first_mean = sum(@view values[1:width]) / width
    last_mean = sum(@view values[(n-width+1):n]) / width
    best_mean = first_mean
    running = sum(@view values[1:width])
    for i = (width+1):n
        running += values[i] - values[i-width]
        best_mean = min(best_mean, max(0.0, running / width))
    end
    # Scales smaller than the useful band do not create unbounded ratios at zero.
    noise_floor = policy.promising_threshold * eps(Float64)
    reduction = 1 - last_mean / max(first_mean, noise_floor)
    regression = last_mean / max(best_mean, noise_floor)
    window = @view values[(n-width+1):n]
    oscillation = maximum(window) / max(minimum(window), policy.promising_threshold)
    enough = n >= 2policy.window
    small = maximum(window) <= policy.promising_threshold
    reduced = reduction >= policy.minimum_reduction
    bounded = regression <= policy.max_regression && oscillation <= policy.max_oscillation
    # A complete window inside the declared useful band needs no further
    # relative reduction. In particular |lambda-1| can cross zero temporarily;
    # comparing a later, still-small window with that accidental historical
    # minimum must not veto the next research Poisson step. Above the band,
    # reduction, regression and oscillation safeguards still all apply.
    # This decision never changes the separate strict acceptance tolerances.
    promising = enough && (small || (reduced && bounded))
    return ResidualTrend(
        Float64(first(values)),
        Float64(last(values)),
        first_mean,
        last_mean,
        best_mean,
        Float64(minimum(values)),
        reduction,
        regression,
        oscillation,
        n > 1 ? Float64(values[end]-values[end-1]) : 0.0,
        promising,
    )
end

"""Assess a history from one fixed Hartree state, using only unmixed map defects.

The caller supplies mathematical usability after solving Dyson and checking
state axes. PSD/current/spectral sum remain independent scientific diagnostics.
A decision before the budget is always :iterate unless the state is unusable.
"""
function research_continuation_decision(
    history::AbstractVector{SCBAIteration},
    options::SolverOptions;
    usable::Bool = true,
)
    trend = options.convergence.trend
    k = residual_trend([row.r_K for row in history], trend)
    sigma = residual_trend([row.r_Σ for row in history], trend)
    lambda = residual_trend([row.r_λ for row in history], trend)
    iteration = isempty(history) ? 0 : last(history).ν
    finite = all(
        row ->
            all(isfinite, (row.r_D, row.r_A, row.r_K, row.r_Σ, row.r_λ, row.λ, row.J)) && row.λ > 0,
        history,
    )
    valid = usable && finite && !isempty(history)
    enabled = options.convergence.mode === :research_continue
    action, reason = if !valid
        (:stop, :unusable_state)
    elseif !enabled
        (:stop, :research_policy_disabled)
    elseif iteration < options.max_scba
        (:iterate, :inner_budget_remaining)
    elseif k.promising && sigma.promising && lambda.promising
        (:continue_poisson, :promising_fixed_point)
    else
        (:stop, :no_convergence_trend)
    end
    return ResearchContinuationDecision(
        action,
        reason,
        iteration,
        options.max_scba,
        valid,
        :last,
        k,
        sigma,
        lambda,
    )
end

function _research_transition_warning(
    scba::SCBAResult,
    options::SolverOptions,
    outer_iteration::Integer,
)
    d = research_continuation_decision(
        scba.history,
        options;
        usable = scba.quality !== :invalid,
    )
    trend = options.convergence.trend
    asdict(x) = Dict(String(name)=>getfield(x, name) for name in fieldnames(typeof(x)))
    return Dict{String,Any}(
        "code"=>d.action === :continue_poisson ? "SCBA_RESEARCH_CONTINUE" :
                "SCBA_RESEARCH_BOUNDED",
        "scope"=>"inner_scba",
        "policy"=>String(options.convergence.mode),
        "trend_contract_version"=>2,
        "outer_iteration"=>Int(outer_iteration),
        "iteration"=>d.iteration,
        "budget"=>d.budget,
        "action"=>String(d.action),
        "reason"=>String(d.reason),
        "selected_state"=>String(d.selected_state),
        "message"=>"Bounded research decision; strict scientific acceptance is independent.",
        "parameters"=>asdict(trend),
        "metrics"=>Dict(
            "r_K"=>asdict(d.keldysh),
            "r_Sigma"=>asdict(d.self_energy),
            "r_lambda"=>asdict(d.normalization),
        ),
    )
end

"""One terminal decision for a fixed-Hartree iteration; no numerical kernels."""
struct SCBAStopDecision
    terminal::Bool
    converged::Bool
    status::Symbol
    quality::Symbol
    reason::String
end

function scba_iteration_decision(
    history::AbstractVector{SCBAIteration},
    options::SolverOptions,
    strict_streak::Integer,
)
    if strict_streak >= options.convergence.required_consecutive_scba_passes
        return SCBAStopDecision(
            true,
            true,
            :converged,
            :strictly_converged,
            "all declared fixed-Hartree strict gates passed",
        )
    end
    if scba_approximate_acceptance(history, options)
        return SCBAStopDecision(
            true,
            false,
            :approximate,
            :approximate_fixed_point,
            "explicit adaptive working band permits Poisson; no strict certificate",
        )
    end
    stopped = _scba_stopping_reason(history, options.tolerances, options.convergence)
    if stopped !== nothing
        return SCBAStopDecision(
            true,
            false,
            stopped,
            :unresolved,
            "strict fixed point stopped before scientific acceptance",
        )
    end
    if !isempty(history) && last(history).ν >= options.max_scba
        continuation = research_continuation_decision(history, options)
        status =
            continuation.action === :continue_poisson ? :research_continue : :max_iterations
        return SCBAStopDecision(
            true,
            false,
            status,
            :unresolved,
            status === :research_continue ?
            "budget exhausted; promising current state permits Poisson with warning" :
            "bounded SCBA budget exhausted without continuation",
        )
    end
    return SCBAStopDecision(false, false, :running, :unresolved, "SCBA budget remains")
end
