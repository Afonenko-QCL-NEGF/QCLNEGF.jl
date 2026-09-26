"""
    ProductionFFTHilbertPlan(ε, Ncolumns;
                              column_chunk=32,
                              worker_count=Threads.nthreads(),
                              verify_roundoff=false)

Reusable, exact-discrete production plan for the zero-padded linear
FFT--Hilbert transform.  Columns are split into fixed, contiguous ranges.
Every logical worker owns a distinct buffer of at most
`L × column_chunk` values and distinct
forward/inverse FFTW plans, so execution has neither shared scratch storage nor
write races.

All FFTW plans are created serially, before any threaded region, with
`num_threads=1`.  FFTW.jl restores its process-wide planner setting after each
plan is made; constructing this object therefore does not change
`FFTW.get_num_threads()`.  The object may be reused by successive SCBA
iterations, but one object must not be used by concurrent outer calls.

Set `verify_roundoff=true` to retain a second set of workers for the doubled
padding length used by `return_roundoff=true`.  This changes only the FFT
summation order and leaves the discrete principal-value map unchanged.

See [the direct Hilbert relation](@ref eq-direct-hilbert),
[Microscopic kernels](@ref theory-kernels), and
[Production backend](@ref theory-production).
"""
struct ProductionFFTHilbertPlan{W}
    energy_count::Int
    column_count::Int
    energy_spacing::Float64
    column_chunk::Int
    primary::W
    verification::Union{Nothing,W}
end

struct _ProductionFFTHilbertWorker{PF,PI}
    columns::UnitRange{Int}
    buffer::Matrix{ComplexF64}
    forward_plan::PF
    inverse_plan::PI
end

struct _ProductionFFTHilbertWorkspace{W}
    transform_length::Int
    kernel_fft::Vector{ComplexF64}
    workers::Vector{W}
end

function _production_fixed_column_ranges(Ncolumns::Integer, worker_count::Integer)
    Ncolumns > 0 || throw(ArgumentError("Ncolumns must be positive"))
    worker_count > 0 || throw(ArgumentError("worker_count must be positive"))
    active_workers = min(Int(Ncolumns), Int(worker_count))
    width, remainder = divrem(Int(Ncolumns), active_workers)
    ranges = Vector{UnitRange{Int}}(undef, active_workers)
    first_column = 1
    for worker = 1:active_workers
        count = width + (worker <= remainder)
        last_column = first_column + count - 1
        ranges[worker] = first_column:last_column
        first_column = last_column + 1
    end
    return ranges
end

function _production_single_thread_plans!(buffer::Matrix{ComplexF64})
    planner_threads = FFTW.get_num_threads()
    forward = FFTW.plan_fft!(buffer, 1; flags = FFTW.ESTIMATE, num_threads = 1)
    inverse = FFTW.plan_ifft!(buffer, 1; flags = FFTW.ESTIMATE, num_threads = 1)
    FFTW.get_num_threads() == planner_threads ||
        error("FFTW planner thread count was not restored")
    return forward, inverse
end

function _production_fft_hilbert_workspace(
    NE::Integer,
    Ncolumns::Integer,
    Δε::Real,
    L::Integer,
    column_chunk::Integer,
    worker_count::Integer,
)
    NE >= 2 || throw(ArgumentError("at least two energy nodes are required"))
    Ncolumns > 0 || throw(ArgumentError("Ncolumns must be positive"))
    Δε > 0 || throw(ArgumentError("energy spacing must be positive"))
    L >= 2NE - 1 || throw(ArgumentError("FFT padding is too short"))
    column_chunk > 0 || throw(ArgumentError("column_chunk must be positive"))
    worker_count > 0 || throw(ArgumentError("worker_count must be positive"))

    planner_threads = FFTW.get_num_threads()
    kernel_fft = zeros(ComplexF64, L)
    for lag = (-(NE-1)):(NE-1)
        iszero(lag) && continue
        kernel_fft[lag+NE] = inv(2π * lag * Δε)
    end
    kernel_plan = FFTW.plan_fft!(kernel_fft, 1; flags = FFTW.ESTIMATE, num_threads = 1)
    mul!(kernel_fft, kernel_plan, kernel_fft)

    active_workers = min(Int(worker_count), Base.Threads.nthreads(:default))
    ranges = _production_fixed_column_ranges(Ncolumns, active_workers)
    buffer_columns = min(Int(column_chunk), maximum(length, ranges))
    workers = map(ranges) do columns
        buffer = zeros(ComplexF64, L, buffer_columns)
        forward, inverse = _production_single_thread_plans!(buffer)
        _ProductionFFTHilbertWorker(columns, buffer, forward, inverse)
    end
    FFTW.get_num_threads() == planner_threads ||
        error("FFTW planner thread count was not restored")
    return _ProductionFFTHilbertWorkspace(Int(L), kernel_fft, workers)
end

function ProductionFFTHilbertPlan(
    ε::AbstractVector{<:Real},
    Ncolumns::Integer;
    column_chunk::Integer = 32,
    worker_count::Integer = Base.Threads.nthreads(:default),
    verify_roundoff::Bool = false,
)
    column_chunk > 0 || throw(ArgumentError("column_chunk must be positive"))
    worker_count > 0 || throw(ArgumentError("worker_count must be positive"))
    NE = length(ε)
    Δε = _uniform_energy_spacing(ε)
    L = nextpow(2, 2NE - 1)
    primary =
        _production_fft_hilbert_workspace(NE, Ncolumns, Δε, L, column_chunk, worker_count)
    verification =
        verify_roundoff ?
        _production_fft_hilbert_workspace(
            NE,
            Ncolumns,
            Δε,
            2L,
            column_chunk,
            worker_count,
        ) : nothing
    return ProductionFFTHilbertPlan(
        NE,
        Int(Ncolumns),
        Δε,
        Int(column_chunk),
        primary,
        verification,
    )
end

function _validate_production_fft_plan(
    plan::ProductionFFTHilbertPlan,
    Γ::AbstractArray{<:Number,4},
    ε::AbstractVector{<:Real},
    wᴱ::AbstractVector{<:Real};
    return_roundoff::Bool = false,
)
    NE, Nk, Nb, Nb2 = size(Γ)
    Nb == Nb2 || throw(DimensionMismatch("Γ blocks must be square"))
    length(ε) == NE == length(wᴱ) || throw(DimensionMismatch("energy data differ"))
    plan.energy_count == NE ||
        throw(DimensionMismatch("FFT plan has the wrong energy size"))
    plan.column_count == Nk * Nb^2 ||
        throw(DimensionMismatch("FFT plan has the wrong column count"))
    Δε = _uniform_energy_spacing(ε)
    isapprox(Δε, plan.energy_spacing; rtol = 256eps(Float64), atol = 0.0) ||
        throw(ArgumentError("FFT plan has a different energy spacing"))
    return_roundoff &&
        isnothing(plan.verification) &&
        throw(ArgumentError("the FFT plan was built without roundoff verification"))
    return nothing
end

function _production_fft_hilbert_apply(
    Γ::AbstractArray{<:Number,4},
    wᴱ::AbstractVector{<:Real},
    workspace::_ProductionFFTHilbertWorkspace,
    column_chunk::Integer,
    task_width::Integer = length(workspace.workers),
)
    task_width >= 1 || throw(ArgumentError("FFT task width must be positive"))
    NE, Nk, Nb, _ = size(Γ)
    Ncolumns = Nk * Nb^2
    Γ2 = reshape(Γ, NE, Ncolumns)
    Λ2 = Matrix{ComplexF64}(undef, NE, Ncolumns)
    L = workspace.transform_length
    kernel_fft = workspace.kernel_fft

    function run_worker(worker_index)
        local worker = workspace.workers[worker_index]
        local buffer = worker.buffer
        for first_column = first(worker.columns):column_chunk:last(worker.columns)
            local last_column = min(first_column + column_chunk - 1, last(worker.columns))
            local active_columns = last_column - first_column + 1
            fill!(buffer, 0.0 + 0.0im)
            @inbounds for local_column = 1:active_columns
                local source_column = first_column + local_column - 1
                for energy = 1:NE
                    buffer[energy, local_column] = wᴱ[energy] * Γ2[energy, source_column]
                end
            end
            mul!(buffer, worker.forward_plan, buffer)
            @inbounds for local_column = 1:active_columns, index = 1:L
                buffer[index, local_column] *= kernel_fft[index]
            end
            mul!(buffer, worker.inverse_plan, buffer)
            @inbounds for local_column = 1:active_columns
                local destination_column = first_column + local_column - 1
                for energy = 1:NE
                    Λ2[energy, destination_column] = buffer[energy+NE-1, local_column]
                end
            end
        end
    end
    active = min(Int(task_width), length(workspace.workers))
    if active == 1
        for worker_index in eachindex(workspace.workers)
            run_worker(worker_index)
        end
    else
        @sync for logical_worker = 1:active
            Base.Threads.@spawn for worker_index =
                    logical_worker:active:length(workspace.workers)

                run_worker(worker_index)
            end
        end
    end
    return reshape(Λ2, size(Γ))
end

"""
    production_fft_hilbert_transform(Γ, ε, wᴱ;
        column_chunk=32, worker_count=Threads.nthreads(),
        return_roundoff=false)
    production_fft_hilbert_transform(Γ, ε, wᴱ, plan;
        return_roundoff=false)

Thread-parallel drop-in production equivalent of
[`fft_hilbert_transform`](@ref).  It evaluates the same zero-padded *linear*
convolution, including the same quadrature weights and omitted principal-value
cell.  Parallelism is only across independent flattened `(k,a,b)` columns;
no physical approximation or reassociation of a reduction is introduced.

The first method constructs plans for a single call.  Production SCBA should
construct one [`ProductionFFTHilbertPlan`](@ref) before its fixed-point loop
and use the second method on every iteration.  A plan is intentionally not
reentrant: concurrent calls require distinct plan objects.

See [the direct Hilbert relation](@ref eq-direct-hilbert),
[the inner SCBA loop](@ref theory-scba), and
[Production backend](@ref theory-production).
"""
function production_fft_hilbert_transform(
    Γ::AbstractArray{<:Number,4},
    ε::AbstractVector{<:Real},
    wᴱ::AbstractVector{<:Real};
    column_chunk::Integer = 32,
    worker_count::Integer = Base.Threads.nthreads(:default),
    return_roundoff::Bool = false,
)
    Ncolumns = size(Γ, 2) * size(Γ, 3) * size(Γ, 4)
    plan = ProductionFFTHilbertPlan(
        ε,
        Ncolumns;
        column_chunk,
        worker_count,
        verify_roundoff = return_roundoff,
    )
    return production_fft_hilbert_transform(Γ, ε, wᴱ, plan; return_roundoff)
end

function production_fft_hilbert_transform(
    Γ::AbstractArray{<:Number,4},
    ε::AbstractVector{<:Real},
    wᴱ::AbstractVector{<:Real},
    plan::ProductionFFTHilbertPlan;
    return_roundoff::Bool = false,
    worker_count::Integer = 0,
)
    _validate_production_fft_plan(plan, Γ, ε, wᴱ; return_roundoff)
    worker_count >= 0 || throw(ArgumentError("FFT worker count cannot be negative"))
    task_width = worker_count == 0 ? length(plan.primary.workers) : Int(worker_count)
    Λ = _production_fft_hilbert_apply(Γ, wᴱ, plan.primary, plan.column_chunk, task_width)
    return_roundoff || return Λ
    Λsecond = _production_fft_hilbert_apply(
        Γ,
        wᴱ,
        something(plan.verification),
        plan.column_chunk,
        task_width,
    )
    # The verification transform is disposable.  Reuse it for the difference
    # so roundoff monitoring does not allocate another full NE×Nk×Nb×Nb field.
    Λsecond .-= Λ
    residual = norm(Λsecond) / (norm(Λ) + 1e-14)
    return Λ, residual
end

"""
    production_fft_workspace_bytes(plan; verification=false)

Exact byte count of explicitly allocated complex worker buffers and the
frequency-domain kernel in `plan`.  FFTW's small opaque plan metadata and the
returned Hilbert array are not included.

See [numerical cost and memory](@ref theory-cost) and
[Production backend](@ref theory-production).
"""
function production_fft_workspace_bytes(
    plan::ProductionFFTHilbertPlan;
    verification::Bool = false,
)
    workspace = verification ? something(plan.verification) : plan.primary
    bytes = sizeof(eltype(workspace.kernel_fft)) * length(workspace.kernel_fft)
    for worker in workspace.workers
        bytes += sizeof(eltype(worker.buffer)) * length(worker.buffer)
    end
    return bytes
end
