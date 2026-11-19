"""
    AlgorithmOptions

Неизменяемые настройки численной реализации и явно выбранных приближений.
Численные ядра используют этот контракт без доступа к реестру, форматам
хранения или управлению запуском.

Первые шесть полей выбирают реализации одних и тех же дискретных уравнений.
При `contraction=:low_rank` и `low_rank_relative_tolerance > 0` применяется
контролируемая численная аппроксимация. Для `kernel_build=:tabulated` измеряется
ошибка интерполяции. Режим `elastic_retarded=:direct_green` задаёт прямое замыкание
упругого SCBA `SigmaR=K[GR]` вместо восстановления через преобразование Гильберта
на конечном окне. Это самостоятельный численный вариант. Последние три переключателя
явно меняют физическую модель SCBA.

См. [выбор метода оптимизации](@ref optimization-decision-tree) и
[проверки физической эквивалентности](@ref theory-validation).
"""
Base.@kwdef struct AlgorithmOptions
    solver_backend::Symbol = :production
    energy_shift::Symbol = :sparse_plan
    hilbert::Symbol = :fft
    contraction::Symbol = :dense_blas
    kernel_build::Symbol = :tabulated
    mixing::Symbol = :linear
    localization::Symbol = :pzp
    embedding::Symbol = :full_cell_resolvent
    embedding_periods::Int = 8
    occupation_normalization::Symbol = :paired_convex
    low_rank_relative_tolerance::Float64 = 0.0
    low_rank_maximum_rank::Int = 0
    anderson_history_depth::Int = 5
    anderson_damping::Float64 = 1.0
    anderson_regularization::Float64 = 1e-12
    elastic_retarded::Symbol = :hilbert
    retarded_real_part::Symbol = :kramers_kronig
    self_energy_structure::Symbol = :full
    transverse_momentum::Symbol = :resolved
end

"""Проверяет настройки численных алгоритмов без обращения к реестру."""
function _check_algorithm_options(options::AlgorithmOptions)
    options.solver_backend in (:educational, :production) ||
        throw(ArgumentError("solver_backend must be :educational or :production"))
    options.energy_shift in (:dense, :sparse_plan, :conservative_pair) || throw(
        ArgumentError("energy_shift must be :dense, :sparse_plan, or :conservative_pair"),
    )
    options.hilbert in (:direct, :fft, :product_integration) ||
        throw(ArgumentError("hilbert must be :direct, :fft, or :product_integration"))
    options.contraction in (:literal, :dense_blas, :low_rank) ||
        throw(ArgumentError("contraction must be :literal, :dense_blas, or :low_rank"))
    options.kernel_build in (:direct, :tabulated, :adaptive_direct) || throw(
        ArgumentError("kernel_build must be :direct, :tabulated, or :adaptive_direct"),
    )
    options.mixing in (:linear, :anderson) ||
        throw(ArgumentError("mixing must be :linear or :anderson"))
    options.localization in (:pzp, :none, :real_space) ||
        throw(ArgumentError("localization must be :pzp, :none, or :real_space"))
    options.embedding in (:full_cell_resolvent, :finite_chain) ||
        throw(ArgumentError("embedding must be :full_cell_resolvent or :finite_chain"))
    options.embedding_periods >= 1 ||
        throw(ArgumentError("embedding_periods must be positive"))
    options.occupation_normalization in (:paired_convex, :scalar_lesser) || throw(
        ArgumentError("occupation_normalization must be :paired_convex or :scalar_lesser"),
    )
    options.retarded_real_part in (:kramers_kronig, :drop) ||
        throw(ArgumentError("retarded_real_part must be :kramers_kronig or :drop"))
    options.elastic_retarded in (:hilbert, :direct_green) ||
        throw(ArgumentError("elastic_retarded must be :hilbert or :direct_green"))
    if options.elastic_retarded === :direct_green &&
       options.retarded_real_part !== :drop && options.contraction === :low_rank &&
       (options.low_rank_relative_tolerance > 0 || options.low_rank_maximum_rank > 0)
        throw(ArgumentError(
            "direct_green elastic retarded closure requires an untruncated kernel; " *
            "general truncated SVD does not establish adjoint/positivity preservation",
        ))
    end
    options.self_energy_structure in (:full, :diagonal) ||
        throw(ArgumentError("self_energy_structure must be :full or :diagonal"))
    options.transverse_momentum in (:resolved, :averaged) ||
        throw(ArgumentError("transverse_momentum must be :resolved or :averaged"))
    0.0 <= options.low_rank_relative_tolerance < 1.0 ||
        throw(ArgumentError("low_rank_relative_tolerance must lie in [0,1)"))
    options.low_rank_maximum_rank >= 0 ||
        throw(ArgumentError("low_rank_maximum_rank cannot be negative"))
    options.anderson_history_depth >= 1 ||
        throw(ArgumentError("anderson_history_depth must be positive"))
    0.0 < options.anderson_damping <= 1.0 ||
        throw(ArgumentError("anderson_damping must lie in (0,1]"))
    options.anderson_regularization >= 0.0 ||
        throw(ArgumentError("anderson_regularization cannot be negative"))
    return options
end


_restart_contract_value(options::ElectronElectronOptions) = Dict(
    String(name)=>_restart_contract_value(value) for
    (name, value) in pairs(electron_electron_identity(options))
)
_restart_contract_value(value::Symbol) = String(value)
_restart_contract_value(value::Union{Real,AbstractString,Nothing}) = value
_restart_contract_value(value) = Dict{String,Any}(
    String(name) => _restart_contract_value(getfield(value, name)) for
    name in fieldnames(typeof(value))
)

"""Численная идентичность для точного продолжения. Пределы итераций задают оставшийся бюджет."""
function _solver_restart_contract(options::SolverOptions, algorithms::AlgorithmOptions)
    return Dict{String,Any}(
        "schema"=>"qcl-negf-exact-restart-v1",
        "algorithms"=>_restart_contract_value(algorithms),
        "solver"=>Dict{String,Any}(
            String(name)=>_restart_contract_value(getfield(options, name)) for
            name in fieldnames(SolverOptions) if name ∉ (:max_scba, :max_poisson)
        ),
    )
end

function _check_exact_restart_contract(
    contract,
    options::SolverOptions,
    algorithms::AlgorithmOptions,
)
    contract === nothing && throw(
        ArgumentError(
            "exact restart is unavailable: checkpoint has no numerical algorithm contract",
        ),
    )
    contract == _solver_restart_contract(options, algorithms) || throw(
        ArgumentError(
            "exact restart numerical contract differs (algorithms, mixing or quality policy); start a new run",
        ),
    )
    return nothing
end

"""Устойчивое описание физической модели с явными размерностями для планирования и точного восстановления."""
function physical_model_identity(models::PhysicalModelOptions)
    result = Dict{String,Any}(
        String(name)=>_restart_contract_value(getfield(models, name)) for
        name in fieldnames(PhysicalModelOptions)
    )
    result["schema"] = "qcl-negf-physical-models-v1"
    return result
end
