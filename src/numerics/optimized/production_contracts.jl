"""
    EnergyShiftPlan

Разреженное непериодическое отображение `X(ε) -> X(ε + δ)` по
кусочно-линейной схеме. Каждая допустимая строка использует не больше
двух исходных узлов. Остальные строки равны нулю.
План эквивалентен [`build_shift_matrix`](@ref). Для его построения и хранения
требуется `O(N_E)`.

См. [периодичность в поле](@ref theory-periodicity)
и [разреженный сдвиг](@ref theory-production).
"""
struct EnergyShiftPlan <: AbstractMatrix{Float64}
    left::Vector{Int}
    right::Vector{Int}
    weight_right::Vector{Float64}
    valid::BitVector
    δ::Float64
end

"""Скалярный показатель в [`SolverEvent`](@ref), независимый от способа записи."""
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

"""Численное событие для функции обратного вызова приложения."""
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

"""Наблюдение принятой точки без доступа к решателю. Профили принадлежат получателю."""
struct ProductionObservation
    stage::Symbol
    status::Symbol
    outer_iteration::Int
    scba_iteration::Int
    current_density_A_per_m2::Float64
    z_nm::Vector{Float64}
    hartree_eV::Vector{Float64}
    density_per_m3::Vector{Float64}
end

function _production_observation(problem, stage, status, outer, scba, current, hartree, density)
    return ProductionObservation(
        stage,
        status,
        outer,
        isempty(scba.history) ? 0 : last(scba.history).ν,
        current,
        problem.grids.x .* (problem.scales.L₀_m * 1e9),
        hartree .* problem.scales.E₀_eV,
        density ./ problem.scales.L₀_m^3,
    )
end

ProductionObservation(solution::NEGFSolution) = _production_observation(
    solution.problem,
    :terminal,
    solution.status,
    length(solution.outer_history),
    solution.scba,
    Float64(ustrip(u"A/m^2", solution.observables[:electron_flow_current])),
    solution.Uᴴ,
    solution.n,
)

_solver_metrics(::Nothing) = SolverMetric[]
_solver_metrics(metric::SolverMetric) = SolverMetric[metric]
_solver_metrics(metrics) = SolverMetric[metric for metric in metrics]

const CHECKPOINT_NONE = :none
const CHECKPOINT_STATE_ONLY = :state_only
const CHECKPOINT_STATE_AND_MARKERS = :state_and_markers
const CHECKPOINT_STATE_AND_ANALYSIS = :state_and_analysis

"""
    ProductionOptions

Настройки ограниченных по памяти операторов, периодичности контрольных
точек и событий прогресса. Они не меняют физическую модель.
`parallel_backend=:threads` распределяет независимые участки энергии
и матричные блоки по потокам Julia. Для него задайте один поток BLAS.
`:blas` использует одного внешнего исполнителя для свёрток, Дайсона,
Келдыша, граничных собственных энергий и невязок, а операции BLAS/LAPACK
могут использовать свои потоки. БПФ с отдельными буферами и простые
поэлементные вычисления всё ещё могут выполняться в потоках Julia.
`worker_count=0` позволяет планировщику использовать все доступные потоки
Julia, положительное значение ограничивает число логических исполнителей.
`residual_chunk` задаёт размер блока независимых задач по невязкам,
а `phase_timing=true` включает хронометраж этапов каждой итерации.
`capture_psd_witness=false` сохраняет скалярные измерения положительности
и отключает сбор матриц и собственного вектора для свидетеля.

`checkpoint_request` задаёт решение о сохранении при наличии приёмника.
Функция получает `(; stage, iteration, outer_iteration)` и возвращает один
из следующих символов: `:none` пропускает сохранение, `:state_only` сохраняет
согласованное состояние с обычной периодичностью физических маркеров,
`:state_and_markers` сохраняет состояние и запрашивает маркеры SCBA на этой карте.
`:state_and_analysis` запрашивает сохранение и аналитический артефакт на внешней границе.
Соответствующие константы: `QCLNumerics.CHECKPOINT_NONE`,
`QCLNumerics.CHECKPOINT_STATE_ONLY`, `QCLNumerics.CHECKPOINT_STATE_AND_MARKERS`
и `QCLNumerics.CHECKPOINT_STATE_AND_ANALYSIS`.
Без отдельной функции заданная периодичность контрольных точек
запрашивает состояние вместе с маркерами.

См. [оптимизированный решатель](@ref theory-production),
[численную стоимость и память](@ref theory-cost)
и [контракт сохранения состояния](https://github.com/Afonenko-QCL-NEGF/qcl-negf/blob/main/runner/docs/src/user/results.md).
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
    capture_psd_witness::Bool = true
end

ProductionOptions(args::Vararg{Any,18}) = ProductionOptions(args..., true)

"""
Вернуть копию `options` с явно указанными изменениями полей.
Единый способ обновления настроек для конфигураций, серий расчётов,
автоподбора и событий прогресса.
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
    capture_psd_witness = options.capture_psd_witness,
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
        capture_psd_witness,
    )
end

"""
Консервативная оценка постоянной и пиковой памяти и объёма вычислений.

См. [численные масштабы и стоимость](@ref theory-cost)
и [ограничение памяти](@ref theory-production).
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
Общий интерфейс операторов рассеяния для оптимизированного расчёта.

Реализации представляют шестииндексную свёртку напрямую, плотной матрицей,
точным множителем при независимости от импульса либо контролируемым
усечённым SVD. Выбранное представление и его научное влияние фиксирует
[`algorithm_manifest`](@ref QCLNEGF.algorithm_manifest).

См. [оптимизированный решатель](@ref theory-production),
[микроскопические ядра](@ref theory-kernels)
и [выбор оптимизации](@ref optimization-decision-tree).
"""
abstract type AbstractProductionKernel end

"""
Маркер прямой шестииндексной свёртки. [`ProductionCache`](@ref) хранит
исходный нормированный тензор, и прямые циклы используют его без копии.
"""
struct LiteralProductionKernel <: AbstractProductionKernel
    N_k::Int
    N_b::Int
end

"""
Точное плотное отображение рассеяния из `(k′,c,d)` в `(k,a,b)`.

См. [микроскопические ядра](@ref theory-kernels)
и [свёртку BLAS](@ref theory-production).
"""
struct DenseProductionKernel <: AbstractProductionKernel
    matrix::Matrix{ComplexF64}
    N_k::Int
    N_b::Int
end

# Exact rank-one interface vertices and their common finite-angle radial map.
# Constructed only after matching the canonical IFR source recipe; arbitrary
# six-index kernels retain the dense implementation.
struct _SeparableIFRProductionKernel <: AbstractProductionKernel
    left::Matrix{ComplexF64}
    right::Matrix{ComplexF64}
    radial::Matrix{ComplexF64}
    N_k::Int
    N_b::Int
end

"""
Точная сокращённая свёртка для ядра, независимого от обоих радиальных импульсов.

См. [микроскопические ядра](@ref theory-kernels)
и [ветвь с независимым от импульса ядром](@ref theory-production).
"""
struct MomentumIndependentProductionKernel <: AbstractProductionKernel
    block::Matrix{ComplexF64}
    N_k::Int
    N_b::Int
end

"""
Контролируемое представление шестииндексного отображения рассеяния
через усечённое SVD: ``K \\approx U_r S_r V_r^\\dagger``.
`relative_frobenius_residual` — точная норма отбрасываемой части
дискретной матрицы по её сингулярным значениям.

См. [контролируемые численные приближения](@ref optimization-controlled).
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
Неизменяемые планы сдвига и точные операторы ядер для всех итераций.

См. [оптимизированный решатель](@ref theory-production),
[периодичность в поле](@ref theory-periodicity)
и [микроскопические ядра](@ref theory-kernels).
"""
struct ProductionCache{P,H}
    W̃₊ᴱᵖ::EnergyShiftPlan
    W̃₋ᴱᵖ::EnergyShiftPlan
    W̃₊ᴸᴼ::EnergyShiftPlan
    W̃₋ᴸᴼ::EnergyShiftPlan
    kernels::Dict{Symbol,AbstractProductionKernel}
    fft_hilbert_plan::P
    product_hilbert_operator::H
    estimate::ProductionMemoryEstimate
    source_grids::ModelGrids
    source_kernel_arrays::Dict{Symbol,Array{ComplexF64,6}}
    algorithms::AlgorithmOptions
end

"""
Скалярная запись одной рабочей точки. `status` содержит статус внешнего
решения, а `scba_quality` — независимую оценку SCBA: строгая сходимость,
допустимое приближённое решение, отсутствие решения или недопустимость.

См. [серии расчётов и возобновление](@ref theory-production).
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
Результат серии расчётов с краткими сведениями по точкам. Полные состояния
всех точек одновременно в памяти не хранятся.

См. [серии расчётов и возобновление](@ref theory-production)
и [контракт сохранения состояния](https://github.com/Afonenko-QCL-NEGF/qcl-negf/blob/main/runner/docs/src/user/results.md).
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

"""
Запросить решение приложения на безопасной границе до выделения снимка.

В запрос входят только счётчики. Запись и управление остаются в приложении.
Явный запрос задаёт решение на каждой границе. Если функции запроса нет,
действует заданная в численных настройках периодичность.
"""
function _production_checkpoint_request(
    options::ProductionOptions,
    stage::Symbol,
    iteration::Int;
    outer_iteration::Int = options.outer_iteration,
)
    request = options.checkpoint_request
    if request !== nothing
        requested = request((; stage, iteration, outer_iteration))
        requested isa Symbol && (
            requested in (CHECKPOINT_NONE, CHECKPOINT_STATE_ONLY, CHECKPOINT_STATE_AND_MARKERS) ||
            (stage === :outer && requested === CHECKPOINT_STATE_AND_ANALYSIS)
        ) ||
            throw(ArgumentError(
                "unsupported checkpoint request for stage $(stage)",
            ))
        return requested
    end
    cadence =
        stage === :scba ? options.checkpoint_every_scba : options.checkpoint_every_outer
    return cadence > 0 && iteration % cadence == 0 ? CHECKPOINT_STATE_AND_MARKERS :
           CHECKPOINT_NONE
end
