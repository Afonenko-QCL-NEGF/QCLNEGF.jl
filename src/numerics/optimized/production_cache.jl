function _begin_solver_stage(
    options::ProductionOptions,
    stage::Symbol;
    label::AbstractString = "",
    iteration = 0,
    total = nothing,
    metrics = nothing,
)
    sink = options.event_sink
    sink === nothing && return false
    values = _solver_metrics(metrics)
    return sink(
        SolverEvent(:begin, stage, String(label), :running, iteration, total, values, ""),
    ) === true
end

function _end_solver_stage(
    options::ProductionOptions,
    stage::Symbol,
    owned::Bool;
    status::Symbol = :completed,
    iteration = nothing,
    total = nothing,
    metrics = nothing,
    message::AbstractString = "",
)
    owned || return nothing
    sink = options.event_sink
    sink === nothing && return nothing
    values = _solver_metrics(metrics)
    return sink(
        SolverEvent(:end, stage, "", status, iteration, total, values, String(message)),
    )
end

function _update_solver_stage(
    options::ProductionOptions,
    stage::Symbol;
    iteration = nothing,
    total = nothing,
    metrics = SolverMetric[],
    message::AbstractString = "",
)
    sink = options.event_sink
    sink === nothing && return false
    values = _solver_metrics(metrics)
    return sink(
        SolverEvent(
            :progress,
            stage,
            "",
            :running,
            iteration,
            total,
            values,
            String(message),
        ),
    ) === true
end

_solver_current_density(value::Real) = Float64(value) / 1.0e4

function _production_worker_count(backend::Symbol, requested::Integer, jobs::Integer)
    requested ≥ 0 || throw(ArgumentError("worker_count cannot be negative"))
    jobs ≥ 0 || throw(ArgumentError("job count cannot be negative"))
    backend in (:threads, :blas) ||
        throw(ArgumentError("parallel backend must be :threads or :blas"))
    jobs == 0 && return 1
    backend === :blas && return 1
    available = Base.Threads.nthreads(:default)
    wanted = requested == 0 ? available : min(Int(requested), available)
    return max(1, min(wanted, Int(jobs)))
end

_production_worker_count(options::ProductionOptions, jobs::Integer) =
    _production_worker_count(options.parallel_backend, options.worker_count, jobs)

_production_thread_worker_count(options::ProductionOptions, jobs::Integer) =
    _production_worker_count(:threads, options.worker_count, jobs)

function _production_worker_ranges!(count::Integer, workers::Integer, operation::Function)
    count ≥ 0 || throw(ArgumentError("linear work size cannot be negative"))
    count == 0 && return nothing
    function run_worker(worker::Int)
        first_index = fld((worker - 1) * count, workers) + 1
        last_index = fld(worker * count, workers)
        first_index > last_index && return nothing
        operation(worker, first_index, last_index)
        return nothing
    end
    if workers == 1
        run_worker(1)
    else
        Base.Threads.@threads :static for worker = 1:workers
            run_worker(worker)
        end
    end
    return nothing
end

_production_worker_ranges!(operation::Function, count::Integer, workers::Integer) =
    _production_worker_ranges!(count, workers, operation)

function _production_parallel_linear!(count::Integer, options::ProductionOptions, operation)
    workers = _production_thread_worker_count(options, count)
    return _production_worker_ranges!(count, workers) do _, first_index, last_index
        @inbounds for index = first_index:last_index
            operation(index)
        end
    end
end

# Julia lowers `do` syntax by passing the anonymous function as the first
# positional argument.  Keep the count-first method above for ordinary calls
# and provide this order for allocation-free `... do index` call sites.
_production_parallel_linear!(operation, count::Integer, options::ProductionOptions) =
    _production_parallel_linear!(count, options, operation)

Base.size(plan::EnergyShiftPlan) = (length(plan.valid), length(plan.valid))

function Base.getindex(plan::EnergyShiftPlan, i::Int, j::Int)
    @boundscheck checkbounds(plan, i, j)
    plan.valid[i] || return 0.0
    left, right, weight = plan.left[i], plan.right[i], plan.weight_right[i]
    return (j == left ? 1 - weight : 0.0) + (j == right ? weight : 0.0)
end

Base.copy(plan::EnergyShiftPlan) = EnergyShiftPlan(
    copy(plan.left),
    copy(plan.right),
    copy(plan.weight_right),
    copy(plan.valid),
    plan.δ,
)

function Base.:*(plan::EnergyShiftPlan, values::AbstractVector{<:Number})
    length(values) == size(plan, 2) || throw(DimensionMismatch("energy vector differs"))
    result = zeros(promote_type(Float64, eltype(values)), length(values))
    for i in eachindex(plan.valid)
        plan.valid[i] || continue
        left, right, weight = plan.left[i], plan.right[i], plan.weight_right[i]
        result[i] =
            left == right || iszero(weight) ? values[left] :
            (1 - weight) * values[left] + weight * values[right]
    end
    return result
end

# Static validation visits only stored entries. Generic AbstractMatrix reductions
# would scan N_E^2 zeros and undo the purpose of the compact representation.
function _valid_shift_plan(plan::EnergyShiftPlan, N::Int)
    length(plan.left) ==
    length(plan.right) ==
    length(plan.weight_right) ==
    length(plan.valid) ==
    N || return false
    all(isfinite, plan.weight_right) && isfinite(plan.δ) || return false
    return all(
        !plan.valid[i] || (1 <= plan.left[i] <= N && 1 <= plan.right[i] <= N) for i = 1:N
    )
end

function _shift_plan_entries(plan::EnergyShiftPlan, i::Int)
    plan.valid[i] || return ((1, 0.0), (1, 0.0))
    left, right, weight = plan.left[i], plan.right[i], plan.weight_right[i]
    return left == right ? ((left, 1.0), (right, 0.0)) :
           ((left, 1 - weight), (right, weight))
end

function _compact_shift_support_and_mass(plan::EnergyShiftPlan, ε, δ)
    support, mass, nonnegative = 0.0, 0.0, 0.0
    N = length(ε)
    source = 1
    for i in eachindex(ε)
        target = ε[i] + δ
        inside = first(ε) <= target <= last(ε)
        while inside && source < N - 1 && ε[source+1] <= target
            source += 1
        end
        left = inside ? source : 0
        row_mass = 0.0
        for (j, value) in _shift_plan_entries(plan, i)
            row_mass += value
            nonnegative = max(nonnegative, -value)
            if !inside || (j != left && j != left + 1)
                support = max(support, abs(value))
            end
        end
        mass = max(mass, abs(row_mass - (inside ? 1.0 : 0.0)))
    end
    return support, mass, nonnegative
end

function _shift_quadrature_adjoint(plus::EnergyShiftPlan, minus::EnergyShiftPlan, weights)
    defect_squared, scale_squared = 0.0, 0.0
    # Visit the union of the supports of plus and transpose(minus). The second
    # pass includes only entries absent from plus, so no term is counted twice.
    for i in eachindex(weights), (j, value) in _shift_plan_entries(plus, i)
        iszero(value) && continue
        forward, backward = weights[i] * value, weights[j] * minus[j, i]
        defect_squared += abs2(forward - backward)
        scale_squared += abs2(forward) + abs2(backward)
    end
    for j in eachindex(weights), (i, value) in _shift_plan_entries(minus, j)
        iszero(value) && continue
        iszero(plus[i, j]) || continue
        backward = weights[j] * value
        defect_squared += abs2(backward)
        scale_squared += abs2(backward)
    end
    return iszero(scale_squared) ? 0.0 : sqrt(defect_squared / scale_squared)
end

function _shift_pair_validation(
    plus::EnergyShiftPlan,
    minus::EnergyShiftPlan,
    ε,
    weights,
    δ,
    discretization::Symbol,
)
    invalid = (
        wrap = Inf,
        nonnegative = Inf,
        row_mass = Inf,
        quadrature_adjoint = Inf,
        cell_weights = Inf,
    )
    N = length(ε)
    N >= 2 && length(weights) == N && discretization === :nodal_linear || return invalid
    all(isfinite, ε) &&
    all(>(0), diff(ε)) &&
    isfinite(δ) &&
    all(isfinite, weights) &&
    all(>(0), weights) &&
    _valid_shift_plan(plus, N) &&
    _valid_shift_plan(minus, N) || return invalid
    psupport, pmass, pnegative = _compact_shift_support_and_mass(plus, ε, δ)
    msupport, mmass, mnegative = _compact_shift_support_and_mass(minus, ε, -δ)
    return (
        wrap = max(psupport, msupport),
        nonnegative = max(pnegative, mnegative),
        row_mass = max(pmass, mmass),
        quadrature_adjoint = _shift_quadrature_adjoint(plus, minus, weights),
        cell_weights = 0.0,
    )
end

"""
    build_shift_plan(ε, δ)

Build the `O(N_E)` exact sparse equivalent of `build_shift_matrix(ε, δ)`.
The input energy nodes must be strictly increasing.

See [Field-periodic closure](@ref theory-periodicity) and
[Production backend](@ref theory-production).
"""
function build_shift_plan(ε::AbstractVector{<:Real}, δ::Real)
    length(ε) >= 2 || throw(ArgumentError("at least two energy nodes are required"))
    all(isfinite, ε) && all(>(0), diff(ε)) ||
        throw(ArgumentError("energy nodes must be finite and strictly increasing"))
    isfinite(δ) || throw(ArgumentError("energy displacement must be finite"))
    N = length(ε)
    left = zeros(Int, N)
    right = zeros(Int, N)
    weight_right = zeros(Float64, N)
    valid = falses(N)
    j = 1
    for e = 1:N
        target = ε[e] + δ
        if target < ε[1] || target > ε[end]
            continue
        elseif target == ε[end]
            left[e] = N
            right[e] = N
            valid[e] = true
            continue
        end
        while j < N - 1 && ε[j+1] <= target
            j += 1
        end
        left[e] = j
        right[e] = j + 1
        weight_right[e] = Float64((target - ε[j]) / (ε[j+1] - ε[j]))
        valid[e] = true
    end
    return EnergyShiftPlan(left, right, weight_right, valid, Float64(δ))
end

"""
    apply_energy_shift!(destination, plan, X)

Apply a non-periodic two-point energy interpolation to a four-axis
`(N_E,N_k,N_b,N_b)` field. `destination` must not alias `X`.

See [Field-periodic closure](@ref theory-periodicity) and the
[production sparse-shift operator](@ref theory-production).
"""
function apply_energy_shift!(
    destination::AbstractArray{<:Number,4},
    plan::EnergyShiftPlan,
    X::AbstractArray{<:Number,4};
    options::ProductionOptions = ProductionOptions(),
)
    size(destination) == size(X) ||
        throw(DimensionMismatch("source and destination shapes differ"))
    Base.mightalias(destination, X) &&
        throw(ArgumentError("in-place aliasing is not supported"))
    NE = size(X, 1)
    length(plan.valid) == NE ||
        throw(DimensionMismatch("shift plan has the wrong energy size"))
    _production_parallel_linear!(length(destination), options) do index
        e = (index - 1) % NE + 1
        if !plan.valid[e]
            destination[index] = zero(eltype(destination))
        else
            l, r, θ = plan.left[e], plan.right[e], plan.weight_right[e]
            base = index - e
            destination[index] =
                l == r || θ == 0 ? X[base+l] : (1 - θ) * X[base+l] + θ * X[base+r]
        end
    end
    return destination
end

"""
Allocate and apply an [`EnergyShiftPlan`](@ref).

See [Field-periodic closure](@ref theory-periodicity) and
[Production backend](@ref theory-production).
"""
function apply_energy_shift(plan::EnergyShiftPlan, X::AbstractArray{<:Number,4})
    T = promote_type(Float64, eltype(X))
    destination = zeros(T, size(X))
    return apply_energy_shift!(destination, plan, X)
end

function _is_momentum_independent(K::AbstractArray{<:Number,6})
    Nk = size(K, 1)
    reference = view(K, 1, 1, :, :, :, :)
    for m = 1:Nk, mp = 1:Nk
        view(K, m, mp, :, :, :, :) == reference || return false
    end
    return true
end

function _low_rank_production_kernel(
    matrix::Matrix{ComplexF64},
    Nk::Int,
    Nb::Int,
    options::AlgorithmOptions,
)
    decomposition = svd(matrix; full = false)
    singular_values = decomposition.S
    total_norm = norm(singular_values)
    tolerance = options.low_rank_relative_tolerance
    maximum_rank =
        options.low_rank_maximum_rank == 0 ? length(singular_values) :
        min(options.low_rank_maximum_rank, length(singular_values))
    rank = length(singular_values)
    if total_norm > 0 && tolerance > 0
        tail_squared = sum(abs2, singular_values)
        rank = 0
        for candidate in eachindex(singular_values)
            tail_squared -= abs2(singular_values[candidate])
            rank = candidate
            sqrt(max(tail_squared, 0.0)) <= tolerance * total_norm && break
        end
    end
    rank = min(rank, maximum_rank)
    discarded =
        rank == length(singular_values) ? 0.0 : norm(@view singular_values[(rank+1):end])
    residual = total_norm == 0 ? 0.0 : discarded / total_norm
    residual <= tolerance + 64eps(Float64) || throw(
        ArgumentError(
            "low_rank_maximum_rank=$maximum_rank cannot satisfy the requested " *
            "relative tolerance $tolerance (attained $residual)",
        ),
    )
    left =
        Matrix{ComplexF64}(decomposition.U[:, 1:rank] .* transpose(singular_values[1:rank]))
    right_adjoint = Matrix{ComplexF64}(decomposition.Vt[1:rank, :])
    return LowRankProductionKernel(left, right_adjoint, Nk, Nb, rank, residual)
end

function _production_kernel(
    K::AbstractArray{<:Number,6},
    algorithms::AlgorithmOptions = AlgorithmOptions(),
)
    _check_algorithm_options(algorithms)
    Nk, Nk2, Nb, Nb2, Nb3, Nb4 = size(K)
    Nk == Nk2 && Nb == Nb2 == Nb3 == Nb4 ||
        throw(DimensionMismatch("kernel order must be (m,m′,a,c,d,b)"))
    algorithms.contraction === :literal && return LiteralProductionKernel(Nk, Nb)
    if _is_momentum_independent(K)
        block4 = Array(view(K, 1, 1, :, :, :, :))
        # Rows `(a,b)`, columns `(c,d)`, with the first index varying fastest.
        block = Matrix{ComplexF64}(reshape(permutedims(block4, (1, 4, 2, 3)), Nb^2, Nb^2))
        return MomentumIndependentProductionKernel(block, Nk, Nb)
    end
    # Rows `(m,a,b)`, columns `(m′,c,d)`.  This is the literal six-index map.
    matrix = Matrix{ComplexF64}(
        reshape(permutedims(K, (1, 3, 6, 2, 4, 5)), Nk * Nb^2, Nk * Nb^2),
    )
    algorithms.contraction === :low_rank &&
        return _low_rank_production_kernel(matrix, Nk, Nb, algorithms)
    return DenseProductionKernel(matrix, Nk, Nb)
end

_kernel_operator_bytes(::LiteralProductionKernel) = 0

function _kernel_operator_bytes(op::DenseProductionKernel)
    return sizeof(eltype(op.matrix)) * length(op.matrix)
end

function _kernel_operator_bytes(op::MomentumIndependentProductionKernel)
    return sizeof(eltype(op.block)) * length(op.block)
end


function _kernel_operator_bytes(op::LowRankProductionKernel)
    return sizeof(eltype(op.left)) * (length(op.left) + length(op.right_adjoint))
end

function _production_memory_int(value::Integer, label::AbstractString)
    0 <= value <= typemax(Int) ||
        throw(OverflowError("$label exceeds addressable Int-sized storage"))
    return Int(value)
end

function _estimate_production_memory(
    n::NumericalParameters,
    mechanism_count::Integer,
    dense_mechanism_count::Integer,
    options::ProductionOptions;
    lo_kernel::Symbol = :dense,
    compact_shifts::Bool = options.algorithms.energy_shift === :sparse_plan,
    worker_capacity::Union{Nothing,Integer} = nothing,
)
    _check_production_options(options)
    NE, Nk, Nb = n.N_E, n.N_k, n.N_b
    M = Int(mechanism_count)
    Mdense = Int(dense_mechanism_count)
    0 ≤ Mdense ≤ M || throw(ArgumentError("invalid dense mechanism count"))
    lo_kernel in (:dense, :reduced, :absent) ||
        throw(ArgumentError("lo_kernel must be :dense, :reduced, or :absent"))
    lo_kernel === :dense &&
        Mdense == 0 &&
        throw(ArgumentError("a dense LO kernel requires at least one dense mechanism"))
    lo_kernel === :reduced &&
        M == Mdense &&
        throw(ArgumentError("a reduced LO kernel requires at least one reduced mechanism"))
    C = sizeof(ComplexF64)
    R = sizeof(Float64)
    n4 = _production_memory_int(
        big(NE) * big(Nk) * big(Nb)^2,
        "four-axis state element count",
    )
    n6 = _production_memory_int((big(Nk) * big(Nb)^2)^2, "six-axis kernel element count")
    fft_required = _production_memory_int(2big(NE) - 1, "minimum FFT transform length")
    Lfft = try
        nextpow(2, fft_required)
    catch error
        (error isa OverflowError || error isa DomainError) || rethrow()
        throw(OverflowError("FFT transform length exceeds addressable Int-sized storage"))
    end
    Lfft > 0 ||
        throw(OverflowError("FFT transform length exceeds addressable Int-sized storage"))
    chunk = min(options.energy_chunk, NE)
    worker_capacity === nothing ||
        worker_capacity >= 1 ||
        throw(ArgumentError("planning worker capacity must be positive"))
    # Planning often runs in a single-threaded validation process. An explicit
    # target capacity must not be silently capped by that process's thread pool.
    estimate_workers(backend, jobs) =
        worker_capacity === nothing ?
        _production_worker_count(backend, options.worker_count, jobs) :
        backend === :blas ? 1 : max(1, min(Int(worker_capacity), jobs))
    contraction_workers = estimate_workers(options.parallel_backend, cld(NE, chunk))
    hcolumns = min(options.hilbert_columns, Nk * Nb^2)
    fft_workers = estimate_workers(:threads, Nk * Nb^2)
    fft_buffer_columns = min(hcolumns, cld(Nk * Nb^2, fft_workers))
    residual_blocks = _production_memory_int(big(NE) * big(Nk), "residual block count")
    residual_chunks = cld(residual_blocks, options.residual_chunk)
    residual_workers =
        options.parallel_backend === :blas ? 1 : estimate_workers(:threads, residual_chunks)
    residual_worker_bytes =
        big(residual_workers) * (
            5big(Nb)^2 * C +
            big(Nb) * R +
            max(1, 3big(Nb)) * C +
            max(1, 5big(Nb)) * R +
            sizeof(LinearAlgebra.BlasInt)
        )
    breakdown_big = Dict{Symbol,BigInt}(
        :dense_shift_matrices_in_problem => compact_shifts ? big(0) : 4big(NE)^2 * R,
        :compact_shift_operators_in_problem =>
            compact_shifts ? 4big(NE) * (2sizeof(Int) + 2R) : big(0),
        :product_hilbert_operator =>
            options.algorithms.hilbert === :product_integration ? big(NE)^2 * R : big(0),
        :cavity_workspace =>
            options.algorithms.embedding === :finite_chain ?
            3big(n4)*C +
            big(estimate_workers(options.parallel_backend, NE*Nk)) *
            (8big(Nb)^2*C + 64big(Nb)*C + (big(Nb)+1)*sizeof(LinearAlgebra.BlasInt)) +
            # A row vector and at most two source/weight pairs per energy
            # and displacement on the production uniform grid. Include array
            # headers/slack conservatively, separately from Julia GC headroom.
            (2big(options.algorithms.embedding_periods)+1) * big(NE) * 160 +
            2big(Nk)*big(Nb)^2*C +
            2big(NE)*R : big(0),
        :green_state => 4big(n4) * C + 2big(NE) * Nk * R,
        :dyson_inverse_workspace =>
            big(estimate_workers(options.parallel_backend, NE*Nk)) *
            (2big(Nb)^2*C + 64big(Nb)*C + (big(Nb)+1)*sizeof(LinearAlgebra.BlasInt)),
        :stored_scattering => 3big(M) * n4 * C,
        :stored_embedding => 9big(n4) * C,
        # Simultaneously live at the fixed-point residual: all scattering
        # candidates (3M), plus/minus/total embedding candidates (9), the
        # combined candidate (3), and the current total self-energy (3).
        :candidate_generation => (3big(M) + 15) * n4 * C,
        # A sweep caller retains the preceding solution while its copied
        # self-energies and the next Green state are being formed.  Count the
        # complete previous Green + scattering + embedding state even though
        # some arrays are initially shared; this is the safe branch.
        :retained_warm_start => (4 + 3big(M) + 9) * n4 * C,
        :raw_six_axis_kernels => big(M) * n6 * C,
        :production_kernel_cache =>
            (big(Mdense) * n6 + big(M - Mdense) * big(Nb)^4) * C,
        :energy_shift_plans => 4big(NE) * (2sizeof(Int) + 2R),
        :blas_workspace => big(contraction_workers) * 2 * Nk * big(Nb)^2 * chunk * C,
        :fft_workspace =>
            (options.verify_fft_roundoff ? 3 : 1) *
            big(Lfft) *
            (1 + big(fft_workers) * fft_buffer_columns) *
            C,
        :residual_workspace =>
            4big(residual_blocks) * R + 3big(residual_chunks) * R + residual_worker_bytes,
        :scalar_reduction_workspace => big(residual_blocks) * C + big(NE) * R,
        # First/worst/last in current and retained warm state, plus one freshly
        # measured witness before the retention step replaces an old payload.
        :physics_witness_blocks => 7 * (((15+2big(M))*big(Nb)^2+Nb)*C+14R),
        :anderson_history =>
            options.algorithms.mixing === :anderson ?
            2big(options.algorithms.anderson_history_depth) * 3 * (M + 2) * n4 * C :
            big(0),
        :anderson_guard_workspace =>
            options.algorithms.mixing === :anderson ?
            5big(estimate_workers(options.parallel_backend, NE*Nk))*big(Nb)^2*C :
            big(0),
        :anderson_gram =>
            options.algorithms.mixing === :anderson ?
            4big(options.algorithms.anderson_history_depth)^2*R : big(0),
        # The active history and recycled pool jointly contain at most h
        # pairs. Solver resets do not allocate another h-pair generation.
        :candidate_transform_temporaries => 6big(n4)*C,
    )
    breakdown = Dict{Symbol,Int}(
        key => _production_memory_int(value, "production memory component $key") for
        (key, value) in breakdown_big
    )
    resident_keys = (
        :dense_shift_matrices_in_problem,
        :compact_shift_operators_in_problem,
        :product_hilbert_operator,
        :cavity_workspace,
        :green_state,
        :dyson_inverse_workspace,
        :physics_witness_blocks,
        :stored_scattering,
        :stored_embedding,
        :raw_six_axis_kernels,
        :production_kernel_cache,
        :energy_shift_plans,
        :fft_workspace,
        :residual_workspace,
        :anderson_history,
        :anderson_guard_workspace,
        :anderson_gram,
    )
    resident = _production_memory_int(
        sum(big(breakdown[key]) for key in resident_keys),
        "estimated resident memory",
    )
    peak = _production_memory_int(
        big(resident) +
        breakdown[:candidate_generation] +
        breakdown[:candidate_transform_temporaries] +
        breakdown[:retained_warm_start] +
        breakdown[:blas_workspace] +
        breakdown[:scalar_reduction_workspace],
        "estimated peak memory",
    )
    # Every mechanism contains lesser and greater contractions.  LO has two
    # shifted contributions in each component, hence two additional maps.
    # Eight is a conventional real-FLOP equivalent for one complex
    # multiply-add.
    dense_maps = 2Mdense + (lo_kernel === :dense ? 2 : 0)
    reduced_maps = 2(M - Mdense) + (lo_kernel === :reduced ? 2 : 0)
    dense_flops = Float64(big(dense_maps) * NE * n6 * 8)
    reduced_flops = Float64(big(reduced_maps) * NE * (big(Nk) * Nb^2 + big(Nb)^4) * 8)
    return ProductionMemoryEstimate(
        resident,
        peak,
        options.memory_budget_bytes,
        dense_flops + reduced_flops,
        breakdown,
    )
end

"""
    estimate_production_memory(numerical, mechanism_count=4; ...)

Conservative estimate available before an expensive `NEGFProblem` is built.
Set `dense_mechanism_count` to the expected number of momentum-dependent
mechanisms.  For the default model acoustic is reduced exactly, while LO,
impurity, and IFR remain dense.

See [Numerical scales and cost](@ref theory-cost) and the
[production RAM contract](@ref theory-production).
"""
function estimate_production_memory(
    n::NumericalParameters,
    mechanism_count::Integer = 4;
    dense_mechanism_count::Integer = 3,
    lo_kernel::Symbol = (
        mechanism_count == 0 ? :absent : (dense_mechanism_count == 0 ? :reduced : :dense)
    ),
    options::ProductionOptions = ProductionOptions(),
    worker_capacity::Union{Nothing,Integer} = nothing,
)
    return _estimate_production_memory(
        n,
        mechanism_count,
        dense_mechanism_count,
        options;
        lo_kernel,
        worker_capacity,
    )
end

"""
Estimate the actual enabled mechanism mix of a built problem.

See [Numerical scales and cost](@ref theory-cost) and the
[production RAM contract](@ref theory-production).
"""
function estimate_production_memory(
    problem::NEGFProblem;
    options::ProductionOptions = ProductionOptions(),
    worker_capacity::Union{Nothing,Integer} = nothing,
)
    Mdense = count(
        name -> !_is_momentum_independent(problem.kernels.K[name]),
        problem.kernels.enabled,
    )
    lo_kernel = if :LO ∉ problem.kernels.enabled
        :absent
    elseif _is_momentum_independent(problem.kernels.K[:LO])
        :reduced
    else
        :dense
    end
    return _estimate_production_memory(
        problem.numerical,
        length(problem.kernels.enabled),
        Mdense,
        options;
        lo_kernel,
        worker_capacity,
        compact_shifts = problem.W₊ᴱᵖ isa EnergyShiftPlan,
    )
end

"""
    build_production_cache(problem; options=ProductionOptions())

Build exact sparse shifts and BLAS-ready kernel maps.  A conservative peak
estimate is checked before any duplicate flattened kernel storage is
allocated.  This is a RAM guard, not a runtime guarantee: the estimate's
`kernel_flops_per_candidate` must also be inspected for large dense kernels.

See [Production backend](@ref theory-production),
[numerical cost and memory](@ref theory-cost), and
[Microscopic kernels](@ref theory-kernels).
"""
function build_production_cache(
    problem::NEGFProblem;
    options::ProductionOptions = ProductionOptions(),
)
    _check_production_options(options)
    estimate = estimate_production_memory(problem; options)
    estimate.peak_bytes <= options.memory_budget_bytes || throw(
        ArgumentError(
            "estimated peak memory $(estimate.peak_bytes) bytes exceeds the " *
            "resolved execution budget $(options.memory_budget_bytes) bytes",
        ),
    )
    sp = _scaled_physics(problem.physical, problem.scales)
    ε = problem.grids.ε
    kernels = Dict{Symbol,AbstractProductionKernel}()
    for name in problem.kernels.enabled
        kernels[name] = _production_kernel(problem.kernels.K[name], options.algorithms)
    end
    # Report the post-detection estimate (acoustic/alloy are normally reduced).
    Mdense = count(
        name -> !_is_momentum_independent(problem.kernels.K[name]),
        problem.kernels.enabled,
    )
    lo_kernel =
        !haskey(kernels, :LO) ? :absent :
        (_is_momentum_independent(problem.kernels.K[:LO]) ? :reduced : :dense)
    estimate = _estimate_production_memory(
        problem.numerical,
        length(kernels),
        Mdense,
        options;
        lo_kernel,
        compact_shifts = problem.W₊ᴱᵖ isa EnergyShiftPlan,
    )
    fft_workers = _production_thread_worker_count(
        options,
        problem.numerical.N_k * problem.numerical.N_b^2,
    )
    fft_hilbert_plan = ProductionFFTHilbertPlan(
        problem.grids.ε,
        problem.numerical.N_k * problem.numerical.N_b^2;
        column_chunk = options.hilbert_columns,
        worker_count = fft_workers,
        verify_roundoff = options.verify_fft_roundoff,
    )
    return ProductionCache(
        build_shift_plan(ε, +sp.Eᵖ),
        build_shift_plan(ε, -sp.Eᵖ),
        build_shift_plan(ε, +sp.ħωᴸᴼ),
        build_shift_plan(ε, -sp.ħωᴸᴼ),
        kernels,
        fft_hilbert_plan,
        options.algorithms.hilbert === :product_integration ?
        product_integration_hilbert_matrix(problem.grids.ε) : nothing,
        estimate,
        problem.grids,
        problem.kernels.K,
        options.algorithms,
    )
end

function _cache_source_contract(cache::ProductionCache, problem::NEGFProblem)
    cache.source_grids === problem.grids || throw(
        ArgumentError(
            "production cache belongs to a different grid object; rebuild it or " *
            "use retarget_problem/retarget_production_cache",
        ),
    )
    cache.source_kernel_arrays === problem.kernels.K || throw(
        ArgumentError(
            "production cache belongs to different normalized kernels; rebuild it",
        ),
    )
    Set(keys(cache.kernels)) == Set(problem.kernels.enabled) ||
        throw(ArgumentError("production cache mechanisms differ from problem"))
    return nothing
end

function _cache_shift_contract(cache::ProductionCache, problem::NEGFProblem)
    sp = _scaled_physics(problem.physical, problem.scales)
    expected = (+sp.Eᵖ, -sp.Eᵖ, +sp.ħωᴸᴼ, -sp.ħωᴸᴼ)
    plans = (cache.W̃₊ᴱᵖ, cache.W̃₋ᴱᵖ, cache.W̃₊ᴸᴼ, cache.W̃₋ᴸᴼ)
    for (label, plan, δ) in zip(("+Ep", "-Ep", "+LO", "-LO"), plans, expected)
        length(plan.valid) == problem.numerical.N_E ||
            throw(DimensionMismatch("production cache $label plan has wrong length"))
        isapprox(plan.δ, δ; rtol = 256eps(Float64), atol = 0.0) ||
            throw(ArgumentError("production cache $label shift is stale"))
    end
    return nothing
end

function _validated_production_cache(
    cache::ProductionCache,
    problem::NEGFProblem,
    options::ProductionOptions,
)
    selected_discretization =
        options.algorithms.energy_shift === :conservative_pair ?
        :finite_volume_piecewise_constant : :nodal_linear
    problem.energy_shift_discretization === selected_discretization || throw(
        ArgumentError(
            "requested energy-shift algorithm does not match the assembled " *
            "$(problem.energy_shift_discretization) operators; rebuild or retarget the problem " *
            "with the requested energy_shift",
        ),
    )
    _cache_source_contract(cache, problem)
    _cache_shift_contract(cache, problem)
    cache.algorithms.hilbert == options.algorithms.hilbert ||
        throw(ArgumentError("production cache Hilbert mode differs; rebuild the cache"))
    cache.algorithms.contraction == options.algorithms.contraction ||
        throw(ArgumentError("production cache contraction mode differs; rebuild the cache"))
    cache.algorithms.low_rank_relative_tolerance ==
    options.algorithms.low_rank_relative_tolerance || throw(
        ArgumentError("production cache low-rank tolerance differs; rebuild the cache"),
    )
    cache.algorithms.low_rank_maximum_rank == options.algorithms.low_rank_maximum_rank ||
        throw(
            ArgumentError(
                "production cache low-rank rank limit differs; rebuild the cache",
            ),
        )
    Mdense = count(
        name -> !_is_momentum_independent(problem.kernels.K[name]),
        problem.kernels.enabled,
    )
    lo_kernel =
        !haskey(cache.kernels, :LO) ? :absent :
        (_is_momentum_independent(problem.kernels.K[:LO]) ? :reduced : :dense)
    estimate = _estimate_production_memory(
        problem.numerical,
        length(cache.kernels),
        Mdense,
        options;
        lo_kernel,
        compact_shifts = problem.W₊ᴱᵖ isa EnergyShiftPlan,
    )
    plan = cache.fft_hilbert_plan
    expected_columns = problem.numerical.N_k * problem.numerical.N_b^2
    expected_workers = _production_thread_worker_count(options, expected_columns)
    plan.energy_count == problem.numerical.N_E &&
    plan.column_count == expected_columns &&
    plan.column_chunk == options.hilbert_columns &&
    length(plan.primary.workers) == expected_workers &&
    (isnothing(plan.verification) == !options.verify_fft_roundoff) ||
        throw(ArgumentError("production cache FFT options differ; rebuild the cache"))
    estimate.peak_bytes <= options.memory_budget_bytes || throw(
        ArgumentError(
            "estimated peak memory $(estimate.peak_bytes) bytes exceeds the " *
            "resolved execution budget $(options.memory_budget_bytes) bytes",
        ),
    )
    return ProductionCache(
        cache.W̃₊ᴱᵖ,
        cache.W̃₋ᴱᵖ,
        cache.W̃₊ᴸᴼ,
        cache.W̃₋ᴸᴼ,
        cache.kernels,
        cache.fft_hilbert_plan,
        cache.product_hilbert_operator,
        estimate,
        cache.source_grids,
        cache.source_kernel_arrays,
        cache.algorithms,
    )
end
