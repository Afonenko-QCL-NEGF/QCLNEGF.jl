mutable struct _AndersonWorkspace
    states::Vector{Vector{SelfEnergyFamily}}
    residuals::Vector{Vector{SelfEnergyFamily}}
    gram::Matrix{Float64}
    free_states::Vector{Vector{SelfEnergyFamily}}
    free_residuals::Vector{Vector{SelfEnergyFamily}}
end

function _AndersonWorkspace(
    states = Vector{SelfEnergyFamily}[],
    residuals = Vector{SelfEnergyFamily}[],
)
    length(states) == length(residuals) ||
        throw(DimensionMismatch("Anderson history differs"))
    count = length(residuals)
    gram = zeros(Float64, count, count)
    for j = 1:count, i = 1:j
        gram[i, j] = gram[j, i] = _anderson_inner(residuals[i], residuals[j])
    end
    return _AndersonWorkspace(
        states,
        residuals,
        gram,
        Vector{SelfEnergyFamily}[],
        Vector{SelfEnergyFamily}[],
    )
end

struct _AndersonGuardWorkspace
    retarded::Matrix{ComplexF64}
    lesser::Matrix{ComplexF64}
    greater::Matrix{ComplexF64}
    hermitian::Matrix{ComplexF64}
    difference::Matrix{ComplexF64}
end
_AndersonGuardWorkspace(n) =
    _AndersonGuardWorkspace((zeros(ComplexF64, n, n) for _ = 1:5)...)

# Same spectral norms, PSD backward-error budget and Keldysh identity as the
# literal guard. Scratch is worker-owned; no norms are weakened for speed.
function _anderson_block_failure!(w, family, e, m)
    copyto!(w.retarded, view(family.Σᴿ, e, m, :, :))
    all(isfinite, w.retarded) || return 1
    for b in axes(w.retarded, 2), a in axes(w.retarded, 1)
        w.hermitian[a, b] = (w.retarded[a, b]-conj(w.retarded[b, a]))/(2im)
    end
    retarded_norm = opnorm(w.retarded)
    eigmax(Hermitian(w.hermitian)) <=
    64eps(Float64)*max(retarded_norm, floatmin(Float64)) || return 2
    copyto!(w.lesser, view(family.Σˡ, e, m, :, :))
    copyto!(w.greater, view(family.Σᵍ, e, m, :, :))
    all(isfinite, w.lesser) && all(isfinite, w.greater) || return 1
    for (block, factor) in ((w.lesser, -im), (w.greater, im))
        # Preserve the literal spectral-norm argument, including its rounding.
        @. w.difference = factor*block
        budget = psd_error_budget(opnorm(w.difference), size(block, 1))
        for b in axes(block, 2), a in axes(block, 1)
            value = factor*block[a, b]
            conjugate = conj(factor*block[b, a])
            w.difference[a, b] = value-conjugate
            w.hermitian[a, b] = (value+conjugate)/2
        end
        norm(w.difference)/2 <= budget || return 3
        eigmin(Hermitian(w.hermitian)) >= -budget || return 4
    end
    for b in axes(w.retarded, 2), a in axes(w.retarded, 1)
        w.difference[a, b] =
            w.retarded[a, b]-conj(w.retarded[b, a])-(w.greater[a, b]-w.lesser[a, b])
    end
    budget = psd_error_budget(
        max(retarded_norm, opnorm(w.lesser), opnorm(w.greater)),
        size(w.retarded, 1),
    )
    norm(w.difference) <= budget || return 5
    return 0
end

function _anderson_guard(families::Vector{SelfEnergyFamily}, options::ProductionOptions)
    isempty(families) &&
        return (passed = true, first_failure = 0, failure_code = 0, checked_blocks = 0)
    ne, nk, nb, _ = size(first(families).Σᴿ)
    blocks = ne*nk*length(families)
    workers = _production_worker_count(options, blocks)
    failures, codes, checked =
        fill(typemax(Int), workers), zeros(Int, workers), zeros(Int, workers)
    function evaluate(worker)
        w = _AndersonGuardWorkspace(nb)
        for index = worker:workers:blocks
            family = families[(index-1)÷(ne*nk)+1]
            local_index = (index-1) % (ne*nk)
            e, m = local_index ÷ nk+1, local_index % nk+1
            code = _anderson_block_failure!(w, family, e, m)
            checked[worker] += 1
            if code != 0
                failures[worker], codes[worker] = index, code
                break
            end
        end
    end
    if workers == 1
        evaluate(1)
    else
        @sync for worker = 1:workers
            Base.Threads.@spawn evaluate(worker)
        end
    end
    index = argmin(failures)
    passed = failures[index] == typemax(Int)
    return (
        passed = passed,
        first_failure = passed ? 0 : failures[index],
        failure_code = codes[index],
        checked_blocks = sum(checked),
    )
end
_anderson_causal(
    families::Vector{SelfEnergyFamily},
    options::ProductionOptions = ProductionOptions(),
) = _anderson_guard(families, options).passed

function _ordered_iteration_families(scattering, plus, minus, order)
    result = SelfEnergyFamily[scattering[name] for name in order]
    push!(result, plus, minus)
    return result
end

function _anderson_residual(candidate::SelfEnergyFamily, state::SelfEnergyFamily)
    return SelfEnergyFamily(
        candidate.Σᴿ .- state.Σᴿ,
        candidate.Σˡ .- state.Σˡ,
        candidate.Σᵍ .- state.Σᵍ,
    )
end

function _anderson_inner(left::Vector{SelfEnergyFamily}, right::Vector{SelfEnergyFamily})
    length(left) == length(right) || throw(DimensionMismatch("Anderson states differ"))
    value = 0.0
    for (a, b) in zip(left, right)
        value += real(dot(vec(a.Σᴿ), vec(b.Σᴿ)))
        value += real(dot(vec(a.Σˡ), vec(b.Σˡ)))
        value += real(dot(vec(a.Σᵍ), vec(b.Σᵍ)))
    end
    return value
end

function _anderson_coefficients(
    residuals::Vector{Vector{SelfEnergyFamily}},
    regularization::Real,
)
    count = length(residuals)
    count >= 1 || throw(ArgumentError("empty Anderson history"))
    count == 1 && return ones(Float64, 1)
    gram = Matrix{Float64}(undef, count, count)
    for j = 1:count, i = 1:j
        gram[i, j] = gram[j, i] = _anderson_inner(residuals[i], residuals[j])
    end
    return _anderson_coefficients(gram, regularization)
end

function _anderson_coefficients(input_gram::Matrix{Float64}, regularization::Real)
    count = size(input_gram, 1)
    count >= 1 || throw(ArgumentError("empty Anderson history"))
    count == 1 && return ones(Float64, 1)
    gram = copy(input_gram)
    scale = max(maximum(abs, diag(gram)), eps(Float64))
    for i = 1:count
        gram[i, i] += regularization * scale
    end
    # The small symmetric Gram system is rank-revealed before solving;
    # near-dependent histories do not enter an unstable indefinite KKT solve.
    decomposition = eigen(Symmetric(gram / scale))
    cutoff = max(regularization, sqrt(eps(Float64)))
    projected = decomposition.vectors' * ones(Float64, count)
    weights = [
        decomposition.values[i] > cutoff ? projected[i] / decomposition.values[i] : 0.0
        for i = 1:count
    ]
    coefficients = decomposition.vectors * weights
    denominator = sum(coefficients)
    isfinite(denominator) && abs(denominator) > eps(Float64) ||
        throw(ArgumentError("Anderson history is rank deficient"))
    coefficients ./= denominator
    all(isfinite, coefficients) ||
        throw(ArgumentError("Anderson coefficient solve produced non-finite values"))
    return coefficients
end

function _anderson_clear!(workspace::_AndersonWorkspace)
    append!(workspace.free_states, workspace.states)
    append!(workspace.free_residuals, workspace.residuals)
    empty!(workspace.states)
    empty!(workspace.residuals)
    workspace.gram = zeros(Float64, 0, 0)
end

function _anderson_append!(workspace, state, candidate, depth)
    prior_norm = isempty(workspace.residuals) ? nothing : workspace.gram[end, end]
    # Evict BEFORE materialization, then reuse its arrays. The history never
    # creates a transient depth+1 full state/residual pair.
    if length(workspace.states) >= depth
        push!(workspace.free_states, popfirst!(workspace.states))
        push!(workspace.free_residuals, popfirst!(workspace.residuals))
        workspace.gram = workspace.gram[2:end, 2:end]
    end
    saved_state =
        isempty(workspace.free_states) ? [_copy_family(f) for f in state] :
        pop!(workspace.free_states)
    residual =
        isempty(workspace.free_residuals) ?
        [_anderson_residual(candidate[j], state[j]) for j in eachindex(state)] :
        pop!(workspace.free_residuals)
    for j in eachindex(state), component in (:Σᴿ, :Σˡ, :Σᵍ)
        copyto!(getfield(saved_state[j], component), getfield(state[j], component))
        out = getfield(residual[j], component)
        raw, current = getfield(candidate[j], component), getfield(state[j], component)
        @. out = raw-current
    end
    residual_norm = _anderson_inner(residual, residual)
    reset =
        prior_norm !== nothing &&
        residual_norm > (1+64eps(Float64))*max(prior_norm, floatmin(Float64))
    reset && _anderson_clear!(workspace)
    count = length(workspace.states)+1
    gram = zeros(Float64, count, count)
    count > 1 && (gram[1:(end-1), 1:(end-1)] .= workspace.gram)
    for i = 1:(count-1)
        gram[i, count] = gram[count, i] = _anderson_inner(workspace.residuals[i], residual)
    end
    gram[count, count] = residual_norm
    workspace.gram = gram
    push!(workspace.states, saved_state)
    push!(workspace.residuals, residual)
    return reset
end

function _anderson_mix!(
    workspace::_AndersonWorkspace,
    state::Vector{SelfEnergyFamily},
    candidate::Vector{SelfEnergyFamily},
    algorithms::AlgorithmOptions,
    α::Real,
    options::ProductionOptions,
)
    length(state) == length(candidate) ||
        throw(DimensionMismatch("Anderson state and candidate differ"))
    copy_started = time_ns()
    reset =
        _anderson_append!(workspace, state, candidate, algorithms.anderson_history_depth)
    history_seconds = (time_ns()-copy_started)*1e-9
    reset_reason = reset ? "residual_increase" : "none"
    coefficient_started = time_ns()
    coefficients = try
        _anderson_coefficients(workspace.gram, algorithms.anderson_regularization)
    catch error
        error isa InterruptException && rethrow()
        reset_reason = "coefficient_solve"
        [zeros(Float64, length(workspace.residuals)-1); 1.0]
    end
    if sum(abs, coefficients) > 10 || !all(isfinite, coefficients)
        reset_reason = "coefficient_bound"
        coefficients = [zeros(Float64, length(workspace.residuals)-1); 1.0]
    end
    coefficient_seconds = (time_ns()-coefficient_started)*1e-9
    depth = length(workspace.states)
    damping = α*algorithms.anderson_damping
    for family_index in eachindex(state), component in (:Σᴿ, :Σˡ, :Σᵍ)
        output = getfield(state[family_index], component)
        _production_parallel_linear!(length(output), options) do index
            value = zero(ComplexF64)
            for history_index in eachindex(coefficients)
                x =
                    getfield(workspace.states[history_index][family_index], component)[index]
                r =
                    getfield(workspace.residuals[history_index][family_index], component)[index]
                value += coefficients[history_index]*(x+damping*r)
            end
            output[index] = value
        end
    end
    guard_started = time_ns()
    guard = _anderson_guard(state, options)
    accelerated_accepted = guard.passed && count(!iszero, coefficients) > 1
    if !guard.passed
        current, raw = last(workspace.states), last(workspace.residuals)
        for j in eachindex(state), component in (:Σᴿ, :Σˡ, :Σᵍ)
            output = getfield(state[j], component)
            x, r = getfield(current[j], component), getfield(raw[j], component)
            @. output = x+damping*r
        end
        _anderson_clear!(workspace)
        reset_reason = "causality_guard"
        _anderson_causal(state, options) || throw(
            DomainError(
                :invalid_keldysh_candidate,
                "convex Anderson fallback is not admissible; inspect raw self-energies",
            ),
        )
    end
    _update_solver_stage(
        options,
        :scba;
        metrics = SolverMetric[
            SolverMetric(:anderson_history_depth, depth),
            SolverMetric(:anderson_reset_reason, reset_reason),
            SolverMetric(:anderson_coefficients_l1, sum(abs, coefficients)),
            SolverMetric(:anderson_acceleration_accepted, accelerated_accepted),
            SolverMetric(:anderson_guard_checked_blocks, guard.checked_blocks),
            SolverMetric(:anderson_guard_first_failure, guard.first_failure),
            SolverMetric(:anderson_guard_failure_code, guard.failure_code),
            SolverMetric(:t_anderson_history_gram, history_seconds, "s"),
            SolverMetric(:t_anderson_coefficients, coefficient_seconds, "s"),
            SolverMetric(:t_anderson_guard, (time_ns()-guard_started)*1e-9, "s"),
        ],
        message = "Anderson bounded history and exact physical safeguard",
    )
    return state
end
