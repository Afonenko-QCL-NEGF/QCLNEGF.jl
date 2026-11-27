"""
    Layer(; d, material, x_Al, Eᶜ, mᶻᵣ, m_parallelᵣ, εᵣ, doped=false)

Один эпитаксиальный слой. Толщину `d` и энергию `Eᶜ` задают величинами Unitful.
Конструктор переводит их в метры и электронвольты. Относительные массы и
диэлектрическая проницаемость безразмерны.

Физический смысл и порядок слоёв приведены в [паспорте эталонной структуры](@ref theory-reference2019).
Единицы входных величин перечислены в [контракте входных данных](@ref full-input-contract).
"""
struct Layer
    d::LengthQuantity
    material::Symbol
    x_Al::Float64
    Eᶜ::EnergyQuantity
    mᶻᵣ::Float64
    m_parallelᵣ::Float64
    εᵣ::Float64
    doped::Bool
    function Layer(d, material, x_Al, Eᶜ, mᶻᵣ, m_parallelᵣ, εᵣ, doped)
        isfinite(d) && d>0u"m" ||
            throw(ArgumentError("layer thickness must be finite and positive"))
        isfinite(Eᶜ) || throw(ArgumentError("band edge must be finite"))
        isfinite(x_Al) && 0<=x_Al<=1 ||
            throw(ArgumentError("alloy fraction must lie in [0,1]"))
        all(x->isfinite(x)&&x>zero(x), (mᶻᵣ, m_parallelᵣ, εᵣ)) ||
            throw(ArgumentError("masses and permittivity must be finite and positive"))
        new(
            _length(d),
            material,
            Float64(x_Al),
            _energy(Eᶜ),
            Float64(mᶻᵣ),
            Float64(m_parallelᵣ),
            Float64(εᵣ),
            doped,
        )
    end

end

function Layer(;
    d,
    material::Symbol,
    x_Al::Real,
    Eᶜ,
    mᶻᵣ::Real,
    m_parallelᵣ::Real,
    εᵣ::Real,
    doped::Bool = false,
)
    d > 0u"m" || throw(ArgumentError("layer thickness must be positive"))
    0 ≤ x_Al ≤ 1 || throw(ArgumentError("x_Al must lie in [0,1]"))
    mᶻᵣ > 0 && m_parallelᵣ > 0 || throw(ArgumentError("effective masses must be positive"))
    εᵣ > 0 || throw(ArgumentError("relative permittivity must be positive"))
    return Layer(
        _length(d),
        material,
        Float64(x_Al),
        _energy(Eᶜ),
        Float64(mᶻᵣ),
        Float64(m_parallelᵣ),
        Float64(εᵣ),
        doped,
    )
end

"""
Размерные и безразмерные параметры одной слоистой транспортной модели.

См. [модель и соглашения](@ref theory-model),
[паспорт эталонной структуры](@ref theory-reference2019) и
[полный контракт входных данных](@ref full-input-contract).
"""
struct PhysicalParameters
    layers::Vector{Layer}
    N_dop²ᴰ::SheetDensityQuantity
    f_ion::Float64
    F_bias::FieldQuantity
    z₀::LengthQuantity
    E_ref::EnergyQuantity
    Tᴸ::TemperatureQuantity
    Tᴸᴼ::TemperatureQuantity
    ε_s::Float64
    ε_∞::Float64
    ħωᴸᴼ::EnergyQuantity
    q_s::WaveNumberQuantity
    qᴸᴼ_s::WaveNumberQuantity
    Δᴵᶠᴿ::LengthQuantity
    Λᴵᶠᴿ::LengthQuantity
    Ξ::EnergyQuantity
    ρ_m::MassDensityQuantity
    v_s::SpeedQuantity
    ΔV_alloy::Union{Nothing,EnergyQuantity}
    Ω₀::Union{Nothing,VolumeQuantity}
    g_s::Int
    interfaces::Vector{LengthQuantity}
    function PhysicalParameters(
        layers,
        N_dop²ᴰ,
        f_ion,
        F_bias,
        z₀,
        E_ref,
        Tᴸ,
        Tᴸᴼ,
        ε_s,
        ε_∞,
        ħωᴸᴼ,
        q_s,
        qᴸᴼ_s,
        Δᴵᶠᴿ,
        Λᴵᶠᴿ,
        Ξ,
        ρ_m,
        v_s,
        ΔV_alloy,
        Ω₀,
        g_s,
        interfaces,
    )
        !isempty(layers) || throw(ArgumentError("physical model needs at least one layer"))
        isfinite(N_dop²ᴰ) && N_dop²ᴰ>0u"m^-2" ||
            throw(ArgumentError("doping must be finite and positive"))
        isfinite(f_ion) && 0<f_ion<=1 ||
            throw(ArgumentError("ionization fraction must lie in (0,1]"))
        isfinite(F_bias) && F_bias>=0u"V/m" ||
            throw(ArgumentError("bias field must be finite and nonnegative"))
        all(isfinite, (z₀, E_ref)) ||
            throw(ArgumentError("energy and position reference must be finite"))
        all(
            x->isfinite(x)&&x>zero(x),
            (Tᴸ, Tᴸᴼ, ε_s, ε_∞, ħωᴸᴼ, q_s, qᴸᴼ_s, Λᴵᶠᴿ, Ξ, ρ_m, v_s),
        ) || throw(ArgumentError("physical scales must be finite and positive"))
        isfinite(Δᴵᶠᴿ) && Δᴵᶠᴿ>=0u"m" ||
            throw(ArgumentError("IFR height must be finite and nonnegative"))
        (ΔV_alloy===nothing)==(Ω₀===nothing) ||
            throw(ArgumentError("alloy energy and volume must be supplied together"))
        ΔV_alloy===nothing ||
            (isfinite(ΔV_alloy)&&ΔV_alloy>0u"eV"&&isfinite(Ω₀)&&Ω₀>0u"m^3") ||
            throw(ArgumentError("alloy parameters must be finite and positive"))
        g_s>=1 || throw(ArgumentError("spin degeneracy must be positive"))
        all(isfinite, interfaces) ||
            throw(ArgumentError("interface positions must be finite"))
        new(
            layers,
            N_dop²ᴰ,
            f_ion,
            F_bias,
            z₀,
            E_ref,
            Tᴸ,
            Tᴸᴼ,
            ε_s,
            ε_∞,
            ħωᴸᴼ,
            q_s,
            qᴸᴼ_s,
            Δᴵᶠᴿ,
            Λᴵᶠᴿ,
            Ξ,
            ρ_m,
            v_s,
            ΔV_alloy,
            Ω₀,
            g_s,
            interfaces,
        )
    end

end

"""
Размеры и границы численных сеток, запас доверенного окна и затравочное уширение
для прямой реализации.

См. [индексы, сетки и размеры массивов](@ref array-contracts) и
[полный контракт входных данных](@ref full-input-contract).
"""
struct NumericalParameters
    N_z::Int
    N_b::Int
    P_basis::Int
    E_min::EnergyQuantity
    E_max::EnergyQuantity
    N_E::Int
    M_E::EnergyQuantity
    k_max::WaveNumberQuantity
    N_k::Int
    N_φ::Int
    qz_max::WaveNumberQuantity
    N_qz::Int
    η_seed::EnergyQuantity
    function NumericalParameters(
        N_z,
        N_b,
        P_basis,
        E_min,
        E_max,
        N_E,
        M_E,
        k_max,
        N_k,
        N_φ,
        qz_max,
        N_qz,
        η_seed,
    )
        N_z>=3 && N_b>=1 && P_basis>=1 || throw(
            ArgumentError("spatial and basis dimensions must be positive and resolved"),
        )
        N_E>=3 && N_k>=2 && N_φ>=4 && N_qz>=3 ||
            throw(ArgumentError("quadrature dimensions are underspecified"))
        all(isfinite, (E_min, E_max, M_E, k_max, qz_max, η_seed)) ||
            throw(ArgumentError("numerical domain values must be finite"))
        E_max>E_min && M_E>0u"eV" && 2M_E<E_max-E_min ||
            throw(ArgumentError("energy window and trusted margin are inconsistent"))
        k_max>0u"m^-1" && qz_max>0u"m^-1" && η_seed>0u"eV" ||
            throw(ArgumentError("wavevector cutoffs and seed broadening must be positive"))
        new(
            N_z,
            N_b,
            P_basis,
            E_min,
            E_max,
            N_E,
            M_E,
            k_max,
            N_k,
            N_φ,
            qz_max,
            N_qz,
            η_seed,
        )
    end

end

function NumericalParameters(;
    N_z::Integer,
    N_b::Integer,
    P_basis::Integer,
    E_min,
    E_max,
    N_E::Integer,
    M_E,
    k_max,
    N_k::Integer,
    N_φ::Integer,
    qz_max,
    N_qz::Integer,
    η_seed,
)
    N_z ≥ 3 || throw(ArgumentError("N_z must be at least 3"))
    N_b ≥ 1 || throw(ArgumentError("N_b must be positive"))
    P_basis ≥ 1 || throw(ArgumentError("P_basis must be positive"))
    N_E ≥ 3 && N_k ≥ 2 && N_φ ≥ 4 && N_qz ≥ 3 ||
        throw(ArgumentError("all quadratures are underspecified"))
    E_max > E_min || throw(ArgumentError("E_max must exceed E_min"))
    M_E > 0u"eV" || throw(ArgumentError("M_E must be positive"))
    2M_E < E_max - E_min ||
        throw(ArgumentError("the Hilbert trust window must be nonempty"))
    k_max > 0u"m^-1" && qz_max > 0u"m^-1" ||
        throw(ArgumentError("momentum cutoffs must be positive"))
    η_seed > 0u"eV" || throw(ArgumentError("η_seed must be positive"))
    return NumericalParameters(
        Int(N_z),
        Int(N_b),
        Int(P_basis),
        _energy(E_min),
        _energy(E_max),
        Int(N_E),
        _energy(M_E),
        _wavenumber(k_max),
        Int(N_k),
        Int(N_φ),
        _wavenumber(qz_max),
        Int(N_qz),
        _energy(η_seed),
    )
end

"""
    SolverTolerances(; ...)

Безразмерные пороги приёмки для внутреннего отображения SCBA, внешнего
отображения Пуассона, законов сохранения, хвостов на краях сетки и проверок
сокращения мнимых частей. Должны выполняться все включённые критерии.
Неконечная невязка всегда означает отказ. Каждый порог проверяет конструктор
`SolverOptions` при включении записи в настройки решателя.

Критерии определены в разделе [проверки решения](@ref theory-validation)
и сведены в таблицу [контракта входных данных](@ref full-input-contract).
"""
Base.@kwdef struct SolverTolerances
    r_D::Float64 = 1e-11
    r_A::Float64 = 1e-9
    r_K::Float64 = 1e-8
    r_Σ::Float64 = 1e-8
    r_λ::Float64 = 1e-8
    r_P::Float64 = 1e-11
    r_U::Float64 = 1e-7
    r_n::Float64 = 1e-7
    r_neutral::Float64 = 1e-10
    r_J::Float64 = 1e-6
    r_C::Float64 = 1e-8
    r_power::Float64 = 1e-4
    r_PSD::Float64 = 1e-10
    r_caus::Float64 = 1e-10
    r_sum::Float64 = 1e-3
    r_obs::Float64 = 1e-4
    r_ζ::Float64 = 1e-10
    r_imag::Float64 = 1e-10
    r_tail::Float64 = 1e-6
    r_edge::Float64 = 1e-6
    r_roundoff::Float64 = 1e-10
end

"""
    DiagnosticQualityPolicy(; enabled=false, ...)

Явная приближённая приёмка для исследовательских расчётов. Сначала пробуют
целевые допуски, затем более грубые. Число строгих попыток после прохождения
целевых допусков ограничено. Алгебраические критерии остаются строгими, а
физические используют собственные заявленные допуски. При приближённой приёмке
`converged=true` не выставляется. Запись и содержащую её политику сходимости
проверяет конструктор `SolverOptions`.
"""
Base.@kwdef struct DiagnosticQualityPolicy
    enabled::Bool = false
    keldysh_threshold::Float64 = 1e-3
    self_energy_threshold::Float64 = 1e-3
    normalization_threshold::Float64 = 1e-4
    required_consecutive_passes::Int = 3
    target_fixed_point_threshold::Float64 = 1e-4
    coarse_wait_iterations::Int = 64
    strict_attempt_iterations::Int = 8
    observable_threshold::Float64 = 1e-3
    positivity_threshold::Float64 = 1e-6
    causality_threshold::Float64 = 1e-10
    outer_threshold::Float64 = 1e-4
end

"""Параметры оценки тенденции свежих невязок r_K, r_Σ и r_λ до смешивания.

Абсолютный порог перспективности характеризует численную пригодность.
После достаточного снижения допускается позднее плато. По окну считают
среднее арифметическое. Регрессия сравнивает последнее окно с лучшим полным
окном. В версии 2 относительную регрессию оценивают за пределами рабочей
полосы. Если окно целиком лежит внутри неё, допускается исследовательское
продолжение.
"""
Base.@kwdef struct ResearchTrendPolicy
    window::Int = 32
    minimum_reduction::Float64 = 0.5
    promising_threshold::Float64 = 1e-4
    max_regression::Float64 = 2.0
    max_oscillation::Float64 = 10.0

    function ResearchTrendPolicy(
        window::Integer,
        minimum_reduction::Real,
        promising_threshold::Real,
        max_regression::Real,
        max_oscillation::Real,
    )
        window >= 2 || throw(ArgumentError("trend window must be at least two"))
        isfinite(minimum_reduction) && 0 < minimum_reduction < 1 ||
            throw(ArgumentError("minimum_reduction must lie in (0,1)"))
        isfinite(promising_threshold) && promising_threshold > 0 ||
            throw(ArgumentError("promising_threshold must be finite and positive"))
        isfinite(max_regression) && max_regression >= 1 ||
            throw(ArgumentError("max_regression must be finite and at least one"))
        isfinite(max_oscillation) && max_oscillation >= 1 ||
            throw(ArgumentError("max_oscillation must be finite and at least one"))
        new(
            Int(window),
            Float64(minimum_reduction),
            Float64(promising_threshold),
            Float64(max_regression),
            Float64(max_oscillation),
        )
    end
end

"""
    ConvergencePolicy(; ...)

Правила остановки, общие для учебного и расчётного решателя. Пороги невязок
задаются в [`SolverTolerances`](@ref). Эта политика определяет, когда принять
последовательность проверок порогов и когда остановить стагнирующую
последовательность SCBA. Запись проверяется при построении `SolverOptions`,
включая прямой позиционный вызов конструктора.

При `stagnation_window=0` ранняя остановка из-за стагнации отключена.
Иначе лучшая нормированная оценка неподвижной точки в более свежей половине
окна должна улучшиться относительно предыдущей половины хотя бы на
`stagnation_relative_improvement`. Оценка учитывает невязки Келдыша, сырой
самоэнергии и нормировки, а также изменения тока и заселённостей. Физическая
допустимость проверяется отдельно, чтобы устойчивый дефект положительности
не скрывал улучшение итерации неподвижной точки.
"""
Base.@kwdef struct ConvergencePolicy
    mode::Symbol = :strict_fail_fast
    trend::ResearchTrendPolicy = ResearchTrendPolicy()
    minimum_scba_iterations::Int = 2
    required_consecutive_scba_passes::Int = 3
    minimum_poisson_iterations::Int = 2
    required_consecutive_poisson_passes::Int = 2
    stagnation_window::Int = 200
    stagnation_relative_improvement::Float64 = 0.01
    diagnostic_quality::DiagnosticQualityPolicy = DiagnosticQualityPolicy()
end

"""
    SolverOptions(; α_Σ=1, α_P=1, max_scba=1000, max_poisson=100, ...)

Настройки линейного смешивания, пределов итераций, допусков проверки,
а также явной политики сходимости и хвостовых окон по энергии и импульсу.
Смешивание одинаково применяется к запаздывающей, меньшей и большей
компонентам каждого семейства самоэнергии и отражается в записях истории.

См. [неподвижную точку SCBA](@ref theory-scba),
[внешний цикл Пуассона](@ref theory-outer-loop) и
[контракт входных данных](@ref full-input-contract).
"""
Base.@kwdef struct SolverOptions
    α_Σ::Float64 = 1.0
    α_P::Float64 = 1.0
    max_scba::Int = 1000
    max_poisson::Int = 100
    tolerances::SolverTolerances = SolverTolerances()
    convergence::ConvergencePolicy = ConvergencePolicy()
    energy_tail_window_fraction::Float64 = 0.05
    momentum_tail_window_fraction::Float64 = 0.10
    function SolverOptions(
        α_Σ,
        α_P,
        max_scba,
        max_poisson,
        tolerances,
        convergence,
        energy_tail_window_fraction,
        momentum_tail_window_fraction,
    )
        return _check_options(
            new(
                α_Σ,
                α_P,
                max_scba,
                max_poisson,
                tolerances,
                convergence,
                energy_tail_window_fraction,
                momentum_tail_window_fraction,
            ),
        )
    end

end

"""
Логические переключатели пяти явно реализованных каналов SCBA.

Микроскопические модели описаны в разделе [ядра рассеяния](@ref theory-kernels).
Их вклад собирает [отображение SCBA](@ref theory-scba).
"""
Base.@kwdef struct ScatteringOptions
    LO::Bool = true
    acoustic::Bool = true
    impurity::Bool = true
    IFR::Bool = true
    alloy::Bool = false
end

"""
Базовые единицы и вычисленные из них безразмерные константы связи.

См. [Unitful и безразмерное ядро](@ref theory-scaling).
"""
struct ScaleSystem
    E₀::EnergyQuantity
    L₀::LengthQuantity
    m_ref::MassQuantity
    E₀_eV::Float64
    E₀_J::Float64
    L₀_m::Float64
    λ_P::Float64
    J₀::typeof(1.0u"A/m^2")
    function ScaleSystem(E₀, L₀, m_ref, E₀_eV, E₀_J, L₀_m, λ_P, J₀)
        all(x->isfinite(x)&&x>zero(x), (E₀, L₀, m_ref, E₀_eV, E₀_J, L₀_m, λ_P, J₀)) ||
            throw(ArgumentError("solver scales must be finite and positive"))
        isapprox(_electronvolts(E₀), E₀_eV; rtol = 64eps(Float64)) ||
            throw(ArgumentError("E0 eV representation differs"))
        isapprox(ustrip(u"J", uconvert(u"J", E₀)), E₀_J; rtol = 64eps(Float64)) ||
            throw(ArgumentError("E0 joule representation differs"))
        isapprox(_metres(L₀), L₀_m; rtol = 64eps(Float64)) ||
            throw(ArgumentError("L0 metre representation differs"))
        new(E₀, L₀, m_ref, E₀_eV, E₀_J, L₀_m, λ_P, J₀)
    end

end

"""
Все узлы и веса квадратур и маска доверенного энергетического окна в
безразмерных переменных.

Внутреннее представление массивов без проверки. Построитель сетки формирует
оси и квадратуры вместе. Прямая сборка не подтверждает соответствие
`NumericalParameters`.

Оси, размеры и квадратуры описаны в разделе
[индексы, сетки и размеры массивов](@ref array-contracts).
"""
struct ModelGrids
    x::Vector{Float64}
    wˣ::Vector{Float64}
    ε::Vector{Float64}
    wᴱ::Vector{Float64}
    trusted_energy::BitVector
    κ::Vector{Float64}
    wᵏ::Vector{Float64}
    qᶻ::Vector{Float64}
    wᑫᶻ::Vector{Float64}
    φ::Vector{Float64}
end

"""
Кусочно заданные профили материалов и ионизованных доноров на сетке
с узлами в центрах ячеек.

Внутреннее представление массивов без проверки. Построитель профиля
согласует сетку с материалами. Для входных данных используйте проверяемый
построитель.

См. [профиль эталонной структуры](@ref theory-reference2019) и
[контракт массивов](@ref array-contracts).
"""
struct MaterialProfiles
    Eᶜ::Vector{Float64}
    mᶻᵣ::Vector{Float64}
    m_parallelᵣ::Vector{Float64}
    εᵣ::Vector{Float64}
    Nᴰ::Vector{Float64}
    x_Al::Vector{Float64}
    layer_index::Vector{Int}
end

"""
Фиксированный локальный базис ячейки PzP с ортогонализацией Лёвдина и
спроецированные одночастичные операторы.

Внутреннее численное представление без проверки. Построитель базиса и
статическая проверка устанавливают размеры, локализацию, ортогональность и невязки.

См. [оператор БенДаниэла–Дьюка и фиксированный базис PzP](@ref theory-basis).
"""
struct BasisData
    localization::Symbol
    Φ::Matrix{ComplexF64}
    χ::Matrix{ComplexF64}
    centres::Vector{Float64}
    spreads::Vector{Float64}
    overlap_eigenvalues::Vector{Float64}
    κ_overlap::Float64
    H₀::Matrix{ComplexF64}
    Z::Matrix{ComplexF64}
    M⁻¹::Matrix{ComplexF64}
    T₊::Matrix{ComplexF64}
    T₋::Matrix{ComplexF64}
    window_eigenvalues::Vector{Float64}
    Egrid₀::Float64
    r_orth::Float64
    r_eigen::Float64
    r_translation::Float64
    r_translation_nearest::Float64
    r_translation_two_pairs::Float64
end

"""
Нормированные ядра с шестью осями в порядке `(m,m′,a,c,d,b)`.

Внутреннее представление массивов без проверки. Общие оси сетки и базиса и
нормировку устанавливают построители ядер.

В `K[name]` хранится нормированный тензор `Khat_s = Kbar_s/qK_s`, а в
`qᴷ[name]` — положительный скаляр `qK_s = maximum(abs, Kbar_s)`. При
свёртке множитель восстанавливается точно. Это позволяет использовать
слагаемые порядка единицы без изменения отображения.

См. [микроскопические ядра и свёртки](@ref theory-kernels) и
[контракт массивов](@ref array-contracts).
"""
struct KernelSet
    K::Dict{Symbol,Array{ComplexF64,6}}
    qᴷ::Dict{Symbol,Float64}
    Fᴸᴼ::Array{ComplexF64,3}
    enabled::Vector{Symbol}
end

"""Явно задаваемые физические приближения. Параметры вычислительной реализации
не определяют этот выбор.

Экранирование по плотности и температуре использует параболическую двумерную
сжимаемость для потенциала доноров и однородное трёхмерное экранирование
Дебая для объёмных LO-фононов. Приближение Кейна меняет только локальную
поперечную кинетическую энергию. Кинетическое замыкание однородной популяции
горячих LO-фононов требует заданных времени распада и плотности участвующих мод.
"""
Base.@kwdef struct PhysicalModelOptions
    screening::Symbol = :fixed
    screening_temperature_K::Float64 = 0.0
    screening_mass_ratio::Float64 = 0.0
    dispersion::Symbol = :parabolic
    nonparabolicity_per_eV::Float64 = 0.0
    lo_population::Symbol = :thermal
    lo_fixed_occupation::Float64 = 0.0
    lo_decay_ps::Float64 = 0.0
    lo_mode_density_per_m3::Float64 = 0.0
    electron_electron::ElectronElectronOptions = ElectronElectronOptions()

    function PhysicalModelOptions(
        screening::Symbol,
        temperature::Real,
        mass::Real,
        dispersion::Symbol,
        alpha::Real,
        lo_population::Symbol,
        occupation::Real,
        decay::Real,
        mode_density::Real,
        electron_electron::ElectronElectronOptions = ElectronElectronOptions(),
    )
        screening in (:fixed, :density_temperature) ||
            throw(ArgumentError("unsupported screening closure"))
        dispersion in (:parabolic, :kane_inplane) ||
            throw(ArgumentError("unsupported dispersion closure"))
        lo_population in (:thermal, :fixed_occupation, :rate_balance) ||
            throw(ArgumentError("unsupported LO population closure"))
        all(
            x -> isfinite(x) && x >= 0,
            (temperature, mass, alpha, occupation, decay, mode_density),
        ) ||
            throw(ArgumentError("physical model parameters must be finite and nonnegative"))
        screening === :fixed ||
            (temperature > 0 && mass > 0) ||
            throw(
                ArgumentError(
                    "density_temperature screening requires electron temperature and band mass",
                ),
            )
        dispersion === :parabolic &&
            alpha != 0 &&
            throw(ArgumentError("nonzero nonparabolicity requires kane_inplane"))
        lo_population === :rate_balance &&
            !(decay > 0 && mode_density > 0) &&
            throw(
                ArgumentError(
                    "LO rate_balance requires decay time and participating mode density",
                ),
            )
        new(
            screening,
            Float64(temperature),
            Float64(mass),
            dispersion,
            Float64(alpha),
            lo_population,
            Float64(occupation),
            Float64(decay),
            Float64(mode_density),
            electron_electron,
        )
    end
end

"""
Полная безразмерная задача после проверки единиц Unitful и масштабирования.

Позиционный конструктор хранит собранные данные численных расширений и
требует явного выбора дискретизации энергетического сдвига и физических моделей.
Публичные построители собирают согласованные оси, единицы и базис. Перед
физическим расчётом вручную собранную задачу необходимо статически проверить.

Порядок сборки приведён в [плане реализации](@ref implementation-plan).
Поля и размерности перечислены в [контракте входных данных](@ref full-input-contract).
"""
struct NEGFProblem{S<:AbstractMatrix{Float64}}
    physical::PhysicalParameters
    numerical::NumericalParameters
    scattering::ScatteringOptions
    scales::ScaleSystem
    grids::ModelGrids
    profiles::MaterialProfiles
    basis::BasisData
    kernels::KernelSet
    W₊ᴱᵖ::S
    W₋ᴱᵖ::S
    W₊ᴸᴼ::S
    W₋ᴸᴼ::S
    # Discretization belongs to the assembled operators, independently of the
    # algorithm used to apply them (dense multiplication or a sparse plan).
    energy_shift_discretization::Symbol
    models::PhysicalModelOptions
end

"""Координаты компактной истории решателя. Массивы принадлежат получателю."""
struct SolverHistoryGrid
    energy_eV::Vector{Float64}
    momentum_per_nm::Vector{Float64}
    spatial_nodes::Int
    energy_scale_eV::Float64
    length_scale_m::Float64
end

SolverHistoryGrid(problem::NEGFProblem) = SolverHistoryGrid(
    problem.grids.ε .* problem.scales.E₀_eV .+ _electronvolts(problem.physical.E_ref),
    problem.grids.κ ./ problem.scales.L₀_m .* 1e-9,
    problem.numerical.N_z,
    problem.scales.E₀_eV,
    problem.scales.L₀_m,
)

"""
Безразмерные массивы функций Грина размера `(N_E,N_k,N_b,N_b)`.

Внутреннее представление массивов без проверки. Согласованность осей с задачей
и физическую допустимость состояния проверяют решатель и загрузчик контрольных точек.

См. [функции Грина и плотность](@ref theory-greens) и
[контракт массивов](@ref array-contracts).
"""
struct GreenState
    Gᴿ::Array{ComplexF64,4}
    Gˡ::Array{ComplexF64,4}
    Gᵍ::Array{ComplexF64,4}
    A::Array{ComplexF64,4}
    condition_number::Matrix{Float64}
    dyson_scale::Matrix{Float64}
end

"""
Семейство безразмерных запаздывающей, меньшей и большей компонент самоэнергии.

Внутреннее представление изменяемых массивов без проверки. Форму и физический
смысл компонент определяет отображение SCBA.

См. [микроскопические ядра](@ref theory-kernels) и
[отображение SCBA](@ref theory-scba).
"""
mutable struct SelfEnergyFamily
    Σᴿ::Array{ComplexF64,4}
    Σˡ::Array{ComplexF64,4}
    Σᵍ::Array{ComplexF64,4}
end

"""Измеренное свидетельство положительной полуопределённости в единицах решателя.
Компактные матричные блоки для восстановления могут отсутствовать.

Индексы относятся к соответствующей сетке энергии и импульса. При пустых
матрицах сохранено скалярное свидетельство той же итерации. Записываются
матрицы для первой, худшей и последней точек. Скаляр и состояние более
поздней итерации не используются для восстановления этих матриц.
"""
struct SCBAPhysicsWitness
    matrix_kind::Symbol
    energy_index::Int
    momentum_index::Int
    minimum_eigenvalue::Float64
    block_norm::Float64
    backward_error::Float64
    absolute_defect::Float64
    relative_defect::Float64
    ratio::Float64
    hermiticity_defect::Float64
    matrices::Dict{String,Matrix{ComplexF64}}
    eigenvector::Vector{ComplexF64}
end

Base.isequal(a::SCBAPhysicsWitness, b::SCBAPhysicsWitness) = all(
    field->isequal(getfield(a, field), getfield(b, field)),
    fieldnames(SCBAPhysicsWitness),
)

"""Настройки ограниченной диагностики. Критерии приёмки решателя от них не зависят.

Спектральные выборки — детерминированные квантили положительного квадратурного
веса в серединах интервалов. Отрицательные спектральные веса ниже погрешности
округления опускаются только в этом диагностическом распределении. При проверке
положительной полуопределённости используются исходные матрицы.
"""
struct SCBAPhysicsMarkerPolicy
    cadence::Int
    max_spectral_blocks::Int
    relative_mode_weight_floor::Float64
    function SCBAPhysicsMarkerPolicy(;
        cadence::Integer = 10,
        max_spectral_blocks::Integer = 512,
        relative_mode_weight_floor::Real = 1e-10,
    )
        cadence >= 1 || throw(ArgumentError("physics marker cadence must be positive"))
        max_spectral_blocks >= 1 ||
            throw(ArgumentError("spectral sample budget must be positive"))
        0 <= relative_mode_weight_floor < 1 ||
            throw(ArgumentError("relative mode weight floor must lie in [0,1)"))
        new(Int(cadence), Int(max_spectral_blocks), Float64(relative_mode_weight_floor))
    end
end

"""Компонента сырой невязки неподвижной точки до смешивания с тем же
масштабом нормы Фробениуса, что у rΣ."""
struct SCBAChannelMarker
    channel::Symbol
    component::Symbol
    residual_absolute::Float64
    residual_scale::Float64
    residual_relative::Float64
end

"""Интеграл столкновений на одном фиксированном нормированном состоянии Грина
до смешивания SCBA.

Знаковые интегралы частиц и энергии и абсолютные масштабы входящих и
исходящих вкладов безразмерны. Сохранённые мнимые части выявляют погрешность
сокращения и округления. Свежее и смешанное отображения обозначаются разными
значениями `state_kind`.
"""
struct SCBACollisionMarker
    channel::Symbol
    state_kind::Symbol
    particle_signed::Float64
    particle_absolute::Float64
    particle_imaginary::Float64
    energy_signed::Float64
    energy_absolute::Float64
    energy_imaginary::Float64
end

for T in (SCBAChannelMarker, SCBACollisionMarker)
    @eval Base.isequal(a::$T, b::$T) =
        all(field -> isequal(getfield(a, field), getfield(b, field)), fieldnames($T))
end

"""Компактные диагностические измерения одного состояния Грина до смешивания.
Они не участвуют в приёмке.

Заряды выражены безразмерными поверхностными концентрациями. Поправки и
флуктуационно-диссипативное соотношение (FDT) используют матричные нормы
Фробениуса с квадратурными весами. FDT сравнивает матрицы занятых и свободных
состояний при одном химическом потенциале, подобранном по спектральному
целевому числу частиц.

Γ = i(SigmaGreater-SigmaLesser) имеет размерность энергии и описывает уширение.
Отношения по выборочным уширениям служат для диагностики. Для ширины пика
на полувысоте и сходимости квадратуры нужны отдельные проверки.
"""
struct SCBAPhysicalMarkers
    measured_iteration::Int
    raw_hole_charge::Float64
    represented_capacity::Float64
    occupied_fraction_a::Float64
    empty_fraction_c::Float64
    relative_correction_Gn::Float64
    relative_correction_Gp::Float64
    equilibrium_applicable::Bool
    equilibrium_status::Symbol
    equilibrium_mu_eV::Float64
    fdt_raw::Float64
    fdt_normalized::Float64
    equilibrium_abs_current_A_m2::Float64
    delta_energy_eV::Float64
    sampled_gamma_over_dE_q10::Float64
    sampled_gamma_over_dE_q50::Float64
    sampled_gamma_over_dE_q90::Float64
    sampled_spectral_weight_underresolved::Float64
    linewidth_sampled_blocks::Int
    linewidth_status::Symbol
    linewidth_sampling_method::Symbol
    marker_cadence::Int
    linewidth_max_blocks::Int
    relative_mode_weight_floor::Float64
    fresh_map_status::Symbol
    lo_shift_over_dE::Float64
    field_shift_over_dE::Float64
    lo_boundary_occupied_fraction::Float64
    lo_boundary_spectral_fraction::Float64
    field_boundary_occupied_fraction::Float64
    field_boundary_spectral_fraction::Float64
    seed_mu_eV::Float64
    seed_number_ratio::Float64
    fdt_raw_seed_mu::Float64
    fdt_normalized_seed_mu::Float64
    channels::Vector{SCBAChannelMarker}
    collisions::Vector{SCBACollisionMarker}
end

Base.isequal(a::SCBAPhysicalMarkers, b::SCBAPhysicalMarkers) = all(
    field -> isequal(getfield(a, field), getfield(b, field)),
    fieldnames(SCBAPhysicalMarkers),
)

"""
Полный набор невязок одного согласованного состояния SCBA до смешивания.

См. [итерацию внутреннего цикла SCBA](@ref theory-scba) и
[критерии проверки](@ref theory-validation).
"""
struct SCBAIteration
    ν::Int
    r_D::Float64
    r_A::Float64
    r_K::Float64
    r_Σ::Float64
    r_λ::Float64
    λ::Float64
    r_PSD::Float64
    r_caus::Float64
    r_roundoff::Float64
    r_Jchange::Float64
    r_population::Float64
    J::Float64
    raw_charge::Float64
    target_charge::Float64
    normalized_charge::Float64
    lambda_change::Float64
    witness::Union{Nothing,SCBAPhysicsWitness}
    physical_markers::Union{Nothing,SCBAPhysicalMarkers}

end

"""История, которую ведёт алгоритм в принятой [точке цикла SCBA](@ref theory-scba)
перед следующим смешиванием."""
struct SCBAMixerState
    method::Symbol
    states::Vector{Vector{SelfEnergyFamily}}
    residuals::Vector{Vector{SelfEnergyFamily}}
end
SCBAMixerState() =
    SCBAMixerState(:linear, Vector{SelfEnergyFamily}[], Vector{SelfEnergyFamily}[])

"""[Состояние SCBA](@ref theory-scba) при фиксированном потенциале Хартри,
семейства компонент по механизмам рассеяния, сходимость и история алгоритма."""
struct SCBAResult
    green::GreenState
    scattering::Dict{Symbol,SelfEnergyFamily}
    embedding::SelfEnergyFamily
    embedding_plus::SelfEnergyFamily
    embedding_minus::SelfEnergyFamily
    history::Vector{SCBAIteration}
    converged::Bool
    status::Symbol
    quality::Symbol
    restart_contract::Union{Nothing,Dict{String,Any}}
    mixer_state::SCBAMixerState

    function SCBAResult(
        green::GreenState,
        scattering::Dict{Symbol,SelfEnergyFamily},
        embedding::SelfEnergyFamily,
        embedding_plus::SelfEnergyFamily,
        embedding_minus::SelfEnergyFamily,
        history::Vector{SCBAIteration},
        converged::Bool,
        status::Symbol,
        quality::Symbol,
        restart_contract::Union{Nothing,Dict{String,Any}} = nothing,
        mixer_state::SCBAMixerState = SCBAMixerState(),
    )
        _check_inner_scba_result_classification(converged, status, quality)
        return new(
            green,
            scattering,
            embedding,
            embedding_plus,
            embedding_minus,
            history,
            converged,
            status,
            quality,
            restart_contract,
            mixer_state,
        )
    end
end

"""Свежее сырое отображение самоэнергии, вычисленное на одном опубликованном состоянии Грина."""
struct FixedPointAuditCandidate
    scattering::Dict{Symbol,SelfEnergyFamily}
    embedding::SelfEnergyFamily
    plus::SelfEnergyFamily
    minus::SelfEnergyFamily
    roundoff::Float64
end

"""Диагностика ограничения на число частиц в безразмерных единицах поверхностной плотности решателя."""
struct ChargeConstraintDiagnostics
    raw_charge::Float64
    target_charge::Float64
    normalized_charge::Float64
    lambda::Float64
    lambda_change::Float64
    constraint_residual::Float64
    normalization_residual::Float64
end

"""
Одна принятая оценка отображения неподвижной точки внешнего цикла Пуассона.

См. [периодическое уравнение Пуассона и внешний цикл](@ref theory-outer-loop).
"""
struct OuterIteration
    μ::Int
    r_P::Float64
    r_U::Float64
    r_n::Float64
    r_neutral::Float64
    r_J::Float64
    r_Jchange::Float64
    r_population::Float64
    ζ::Float64
    J::Float64
    populations::Vector{Float64}
    density::Vector{Float64}
end

"""
Именованные показатели проверки и сообщения обо всех непройденных критериях приёмки.

См. [проверку решения](@ref theory-validation).
"""
struct ConvergenceReport
    passed::Bool
    metrics::Dict{Symbol,Float64}
    messages::Vector{String}
end

"""
Полный результат замкнутого решателя.

`Uᴴ` и `n` — безразмерные массивы ``U_H/E_0`` и ``L_0^3 n``.
Физические наблюдаемые хранятся в `observables` как величины Unitful.
Статус `status == :converged` выдаётся, когда найдены обе неподвижные точки
и пройдена итоговая независимая проверка.

См. [замкнутый цикл Пуассона–SCBA](@ref theory-outer-loop),
[наблюдаемые и балансы](@ref theory-observables) и
[итоговую проверку](@ref theory-validation).
"""
struct NEGFSolution
    problem::NEGFProblem
    options::SolverOptions
    Uᴴ::Vector{Float64}
    n::Vector{Float64}
    scba::SCBAResult
    outer_history::Vector{OuterIteration}
    observables::Dict{Symbol,Any}
    report::ConvergenceReport
    converged::Bool
    status::Symbol
end

function Base.show(io::IO, problem::NEGFProblem)
    n = problem.numerical
    print(
        io,
        "NEGFProblem( N_z=$(n.N_z), N_b=$(n.N_b), ",
        "N_E=$(n.N_E), N_k=$(n.N_k), mechanisms=$(problem.kernels.enabled))",
    )
end

function Base.show(io::IO, solution::NEGFSolution)
    print(
        io,
        "NEGFSolution(status=$(solution.status), converged=$(solution.converged), ",
        "outer_iterations=$(length(solution.outer_history)))",
    )
end

"""Явно ограниченное расширение области. При каждом запуске операторы строятся заново."""
Base.@kwdef struct DomainAdaptationPolicy
    mode::Symbol = :none
    maximum_expansions::Int = 2
    maximum_energy_nodes::Int = 100_001
    growth_fraction::Float64 = 0.5
    tail_threshold::Float64 = 1e-6
    function DomainAdaptationPolicy(
        mode::Symbol,
        expansions::Integer,
        nodes::Integer,
        growth::Real,
        tail::Real,
    )
        mode in (:none, :expand_energy_window) ||
            throw(ArgumentError("unknown domain adaptation policy"))
        expansions >= 0 || throw(ArgumentError("maximum_expansions cannot be negative"))
        nodes >= 3 || throw(ArgumentError("maximum_energy_nodes must be at least three"))
        isfinite(growth) && growth > 0 ||
            throw(ArgumentError("growth_fraction must be positive"))
        isfinite(tail) && 0 < tail < 1 ||
            throw(ArgumentError("tail_threshold must lie in (0,1)"))
        new(mode, Int(expansions), Int(nodes), Float64(growth), Float64(tail))
    end
end
