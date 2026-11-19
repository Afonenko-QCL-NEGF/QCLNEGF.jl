"""
    reference_parameters(; f_ion=1, F_bias=nothing, V_period=56u"mV",
                           Tᴸ=200u"K", ΔV_alloy=nothing, Ω₀=nothing)

Опубликованная последовательность слоёв эталонной структуры с явно выбранными
параметрами материалов и рассеяния учебной модели. Все размерные входы —
величины Unitful. Несовместимые единицы вызывают `Unitful.DimensionError`.
По умолчанию первично точное падение `V_period=56 mV` и `F_bias=V_period/Lp`.
Явный `F_bias` имеет приоритет над `V_period`. Параметры сплава `ΔV_alloy` и
`Ω₀` либо оба опускают, либо оба задают с размерностями энергии и объёма соответственно.

См. [паспорт эталонной структуры](@ref theory-reference2019),
[соглашения модели](@ref theory-model) и
[полный контракт входных данных](@ref full-input-contract).
"""
function reference_parameters(;
    f_ion::Real = 1.0,
    F_bias = nothing,
    V_period = 56.0u"mV",
    Tᴸ = 200.0u"K",
    Tᴸᴼ = Tᴸ,
    Δᴵᶠᴿ = 0.10u"nm",
    Λᴵᶠᴿ = 9.0u"nm",
    q_s = 0.20u"nm^-1",
    qᴸᴼ_s = 0.20u"nm^-1",
    ΔV_alloy = nothing,
    Ω₀ = nothing,
)

    0 < f_ion ≤ 1 ||
        throw(ArgumentError("f_ion must lie in (0,1] for the fixed-density solver"))
    Tᴸ > 0u"K" && Tᴸᴼ > 0u"K" ||
        throw(ArgumentError("lattice and LO temperatures must be positive"))
    q_s > 0u"m^-1" && qᴸᴼ_s > 0u"m^-1" ||
        throw(ArgumentError("screening wave numbers must be positive"))
    Δᴵᶠᴿ ≥ 0u"m" && Λᴵᶠᴿ > 0u"m" || throw(
        ArgumentError("IFR height must be nonnegative and correlation length positive"),
    )
    (ΔV_alloy === nothing) == (Ω₀ === nothing) ||
        throw(ArgumentError("ΔV_alloy and Ω₀ must be supplied together"))
    if ΔV_alloy !== nothing
        ΔV_alloy > 0u"eV" || throw(ArgumentError("ΔV_alloy must be positive"))
        Ω₀ > 0u"m^3" || throw(ArgumentError("Ω₀ must be positive"))
    end
    barrier =
        (; x_Al = 0.25, Eᶜ = 0.20775u"eV", mᶻᵣ = 0.08775, m_parallelᵣ = 0.08775, εᵣ = 12.9)
    well = (; x_Al = 0.0, Eᶜ = 0.0u"eV", mᶻᵣ = 0.067, m_parallelᵣ = 0.067, εᵣ = 12.9)
    layers = Layer[
        Layer(d = 3.26u"nm", material = :AlGaAs, doped = false; barrier...),
        Layer(d = 7.99u"nm", material = :GaAs, doped = false; well...),
        Layer(d = 1.90u"nm", material = :AlGaAs, doped = false; barrier...),
        Layer(d = 8.40u"nm", material = :GaAs, doped = false; well...),
        Layer(d = 2.90u"nm", material = :GaAs, doped = true; well...),
        Layer(d = 5.16u"nm", material = :GaAs, doped = false; well...),
    ]
    interfaces = _length.([0.0u"nm", 3.26u"nm", 11.25u"nm", 13.15u"nm"])
    Lp = sum(layer.d for layer in layers)
    field = F_bias === nothing ? uconvert(u"V/m", V_period / Lp) : F_bias
    field ≥ 0u"V/m" || throw(ArgumentError("reference design sign convention requires F_bias ≥ 0"))
    return PhysicalParameters(
        layers,
        _sheetdensity(4.5e14u"m^-2"),
        Float64(f_ion),
        _field(field),
        _length(0.0u"nm"),
        _energy(0.0u"eV"),
        _temperature(Tᴸ),
        _temperature(Tᴸᴼ),
        12.9,
        10.89,
        _energy(36.7u"meV"),
        _wavenumber(q_s),
        _wavenumber(qᴸᴼ_s),
        _length(Δᴵᶠᴿ),
        _length(Λᴵᶠᴿ),
        _energy(7.0u"eV"),
        _massdensity(5317.0u"kg/m^3"),
        _speed(5240.0u"m/s"),
        ΔV_alloy === nothing ? nothing : _energy(ΔV_alloy),
        Ω₀ === nothing ? nothing : _volume(Ω₀),
        2,
        interfaces,
    )
end

"""
Эталонная дискретизация из математической постановки.

См. [индексы, сетки и размеры массивов](@ref array-contracts) и
[контракт численных входных данных](@ref full-input-contract).
"""
baseline_numerics() = NumericalParameters(
    N_z = 198,
    N_b = 5,
    P_basis = 3,
    E_min = -0.15u"eV",
    E_max = 0.45u"eV",
    N_E = 2401,
    M_E = 56.0u"meV",
    k_max = 0.60u"nm^-1",
    N_k = 61,
    N_φ = 64,
    qz_max = 10.0u"nm^-1",
    N_qz = 401,
    η_seed = 0.5u"meV",
)

"""
Расчётная дискретизация эталонной структуры. Сохраняются исходные параметры
пространственной, энергетической и импульсной сеток и размер базиса.
Запас `80 meV` для непериодического сдвига и доверенного окна позволяет
проводить штатный расчёт ВАХ до `64 mV` на период внутри сохранённого
энергетического окна. Запас меняет маски валидации. Узлы энергии сохраняются,
а достаточность сетки проверяется отдельно.

См. [расчётный решатель](@ref theory-production),
[индексы, сетки и размеры массивов](@ref array-contracts) и
[численные масштабы и стоимость](@ref theory-cost).
"""
reference_production_numerics() = NumericalParameters(
    N_z = 198,
    N_b = 5,
    P_basis = 3,
    E_min = -0.15u"eV",
    E_max = 0.45u"eV",
    N_E = 2401,
    M_E = 80.0u"meV",
    k_max = 0.60u"nm^-1",
    N_k = 61,
    N_φ = 64,
    qz_max = 10.0u"nm^-1",
    N_qz = 401,
    η_seed = 0.5u"meV",
)

"""
Малая детерминированная сетка для примеров и автоматических проверок. Уравнения
сохраняются, меняется дискретизация. Количественные выводы об эталонной структуре
требуют полноценной сетки.

См. [индексы, сетки и размеры массивов](@ref array-contracts) и обязательные
[проверки и оценку сходимости по сетке](@ref theory-validation).
"""
tutorial_numerics() = NumericalParameters(
    N_z = 48,
    N_b = 3,
    P_basis = 1,
    E_min = -0.10u"eV",
    E_max = 0.28u"eV",
    N_E = 49,
    M_E = 56.0u"meV",
    k_max = 0.45u"nm^-1",
    N_k = 5,
    N_φ = 8,
    qz_max = 6.0u"nm^-1",
    N_qz = 25,
    η_seed = 2.0u"meV",
)

"""
Итоговые пороги приёмки, прямые обновления и установленные защитные пределы.

См. [проверку решения](@ref theory-validation) и
[контракт входных данных неподвижной точки](@ref full-input-contract).
"""
baseline_options() = SolverOptions()

"""
Начальные настройки с демпфированием и строгими допусками для расчёта
неподвижных точек эталонной структуры. Коэффициенты демпфирования —
численные параметры, которые варьируют при исследовании сходимости.
Допуски совпадают с [`baseline_options`](@ref).

См. [неподвижную точку SCBA](@ref theory-scba),
[внешний цикл Пуассона](@ref theory-outer-loop) и
[расчётный решатель](@ref theory-production).
"""
reference_production_solver_options() = SolverOptions(
    α_Σ = 0.10,
    α_P = 0.15,
    max_scba = 2000,
    max_poisson = 100,
    tolerances = SolverTolerances(),
    convergence = ConvergencePolicy(),
)

"""
Мягкие параметры смешивания и пороги для коротких проверочных запусков алгоритма.

См. [неподвижную точку SCBA](@ref theory-scba) и
[проверку решения](@ref theory-validation). Для количественной приёмки
используются строгие допуски.
"""
tutorial_options() = SolverOptions(
    α_Σ = 0.35,
    α_P = 0.25,
    max_scba = 80,
    max_poisson = 30,
    convergence = ConvergencePolicy(
        required_consecutive_scba_passes = 1,
        required_consecutive_poisson_passes = 1,
        stagnation_window = 0,
        stagnation_relative_improvement = 0.0,
    ),
    tolerances = SolverTolerances(
        r_D = 1e-10,
        r_A = 1e-6,
        r_K = 5e-4,
        r_Σ = 5e-4,
        r_λ = 5e-4,
        r_P = 1e-10,
        r_U = 5e-4,
        r_n = 5e-4,
        r_neutral = 1e-8,
        r_J = 1e-3,
        r_C = 1e-3,
        r_power = 5e-2,
        r_PSD = 1e-7,
        r_caus = 1e-7,
        r_sum = 5e-2,
        r_obs = 5e-3,
        r_ζ = 1e-7,
        r_imag = 1e-6,
        r_tail = 1.0,
        r_edge = 1.0,
        r_roundoff = 1e-6,
    ),
)

"""
По умолчанию включено рассеяние на LO-фононах, акустических фононах,
ионизованных примесях и шероховатости интерфейсов (IFR).
Рассеяние на сплаве выключено.

См. [микроскопические ядра](@ref theory-kernels).
"""
default_scattering() = ScatteringOptions()

period_length(p::PhysicalParameters) = sum(layer.d for layer in p.layers)
ionized_sheet_density(p::PhysicalParameters) = p.f_ion * p.N_dop²ᴰ
