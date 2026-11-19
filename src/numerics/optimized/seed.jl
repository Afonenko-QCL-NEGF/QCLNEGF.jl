@doc raw"""
    _seed_number_weights_production(Gᴿ, grids, g_s; thread_trace=nothing)

Предварительно вычислить точный дискретный вклад каждого узла энергии
в поверхностную концентрацию равновесного начального состояния.
При ``A=i(Gᴿ-Gᴬ)`` возвращаемый вектор равен

```math
q_e = \frac{g_s}{2\pi}\,w^E_e
      \sum_m w^k_m\,\operatorname{tr} A_{em},
```

Функционал числа частиц начального состояния принимает вид
``N(\mu)=\sum_e q_e f(E_e;\mu)``. Это та же конечная квадратура, что в
[`number_functional`](@ref). Множитель Ферми суммируется по энергии
после свёртки по импульсу и следу по базису. Для пробного
химпотенциала не нужно выделять массив меньшей компоненты.

`Gᴿ` — плотный комплексный массив размера `(N_E,N_k,N_b,N_b)`.
Результат — вещественный безразмерный вектор длины `N_E`. Независимые
строки энергии обрабатываются в потоках. Внутренний `thread_trace`
длины `N_E` записывает поток каждой строки для проверки расписания.

См. [функции Грина и плотность](@ref theory-greens),
[нормировку числа частиц](@ref eq-number-normalization)
и [начальное состояние](@ref theory-production).
"""
function _seed_number_weights_production(
    Gᴿ::AbstractArray{<:Number,4},
    grids::ModelGrids,
    g_s::Integer;
    options::ProductionOptions = ProductionOptions(),
    thread_trace::Union{Nothing,AbstractVector{<:Integer}} = nothing,
)
    NE, Nk, Nb, Nb2 = size(Gᴿ)
    Nb == Nb2 || throw(DimensionMismatch("Gᴿ blocks must be square"))
    length(grids.ε) == NE && length(grids.wᴱ) == NE ||
        throw(DimensionMismatch("energy grid differs from Gᴿ"))
    length(grids.κ) == Nk && length(grids.wᵏ) == Nk ||
        throw(DimensionMismatch("momentum grid differs from Gᴿ"))
    g_s > 0 || throw(ArgumentError("g_s must be positive"))
    thread_trace === nothing ||
        length(thread_trace) == NE ||
        throw(DimensionMismatch("thread_trace must have length N_E"))

    q = zeros(Float64, NE)
    prefactor = Float64(g_s) / (2π)
    _production_parallel_linear!(NE, options) do e
        # A diagonal element is evaluated with exactly the defining
        # i(Gᴿ-Gᴬ) expression.  Its imaginary roundoff is identically zero,
        # while taking `real` makes the real-valued storage contract explicit.
        local momentum_sum = 0.0
        @inbounds for m = 1:Nk
            local trace_A = 0.0
            for a = 1:Nb
                local z = Gᴿ[e, m, a, a]
                trace_A += real(im * (z - conj(z)))
            end
            momentum_sum += grids.wᵏ[m] * trace_A
        end
        q[e] = prefactor * grids.wᴱ[e] * momentum_sum
        thread_trace === nothing || (thread_trace[e] = Base.Threads.threadid())
    end
    all(isfinite, q) || throw(DomainError(q, "seed spectral number weights must be finite"))
    return q
end

function _seed_number_from_weights(
    μ::Real,
    ε::AbstractVector{<:Real},
    q::AbstractVector{<:Real},
    kBT::Real,
)
    length(ε) == length(q) ||
        throw(DimensionMismatch("energy nodes and seed weights differ"))
    kBT > 0 || throw(ArgumentError("kBT must be positive"))
    value = 0.0
    @inbounds for e in eachindex(ε, q)
        value += q[e] * _fermi(ε[e], μ, kBT)
    end
    return value
end

@doc raw"""
    _seed_lesser_production!(Gˡ, Gᴿ, ε, μ, kBT;
                             thread_trace=nothing)

Заполнить меньшую компоненту равновесного начального состояния
по ``G^<(E_e,k_m)=i f(E_e;\mu) A(E_e,k_m)`` и
``A=i(G^R-G^A)``. Независимые блоки `(e,m)` распределяются по потокам
Julia. Все массивы безразмерны и имеют форму `(N_E,N_k,N_b,N_b)`.

Необязательный внутренний вектор `thread_trace` длины `N_E*N_k`
сохраняет исполнителя каждого блока для проверки расписания.

См. [функции Грина и плотность](@ref theory-greens)
и [оптимизированный решатель](@ref theory-production).
"""
function _seed_lesser_production!(
    Gˡ::AbstractArray{ComplexF64,4},
    Gᴿ::AbstractArray{<:Number,4},
    ε::AbstractVector{<:Real},
    μ::Real,
    kBT::Real;
    options::ProductionOptions = ProductionOptions(),
    thread_trace::Union{Nothing,AbstractVector{<:Integer}} = nothing,
)
    size(Gˡ) == size(Gᴿ) || throw(DimensionMismatch("Gˡ and Gᴿ arrays differ"))
    NE, Nk, Nb, Nb2 = size(Gᴿ)
    Nb == Nb2 || throw(DimensionMismatch("Gᴿ blocks must be square"))
    length(ε) == NE || throw(DimensionMismatch("energy axis differs"))
    kBT > 0 || throw(ArgumentError("kBT must be positive"))
    thread_trace === nothing ||
        length(thread_trace) == NE * Nk ||
        throw(DimensionMismatch("thread_trace must have length N_E*N_k"))

    _production_parallel_linear!(NE * Nk, options) do linear
        local e = (linear - 1) % NE + 1
        local m = (linear - 1) ÷ NE + 1
        local f = _fermi(ε[e], μ, kBT)
        @inbounds for b = 1:Nb, a = 1:Nb
            # i*A = -(Gᴿ-Gᴬ), evaluated directly without a temporary A.
            Gˡ[e, m, a, b] = -f * (Gᴿ[e, m, a, b] - conj(Gᴿ[e, m, b, a]))
        end
        thread_trace === nothing || (thread_trace[linear] = Base.Threads.threadid())
    end
    return Gˡ
end

"""
    _seed_green_production(problem, h,
                           options=ProductionOptions(); thread_trace=nothing)

Построить равновесное начальное состояние с заданным числом частиц.
Многопоточный оператор Дайсона выполняется один раз. При делении
пополам по химпотенциалу используются `N_E` точных дискретных спектральных
весов. Окончательная меньшая компонента `Gˡ` строится параллельно
по блокам `(E_e,k_m)`.

Результат `(green, μ₀)` соответствует контракту `_seed_green`.
У [`GreenState`](@ref) четыре комплексных поля размера `(N_E,N_k,N_b,N_b)`.
`μ₀` — безразмерная энергия. Дискретная квадратура сохраняется.

Для тестов можно передать `thread_trace=(weights=..., lesser=...)`
с целочисленными векторами длины `N_E` и `N_E*N_k`.

См. [неподвижную точку SCBA](@ref theory-scba),
[функции Грина и плотность](@ref theory-greens)
и [оптимизированный решатель](@ref theory-production).
"""
function _seed_green_production(
    problem::NEGFProblem,
    h::AbstractArray{<:Number,3},
    options::ProductionOptions = ProductionOptions();
    thread_trace = nothing,
)
    n = problem.numerical
    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    size(h) == (n.N_k, n.N_b, n.N_b) || throw(DimensionMismatch("h has wrong shape"))
    _check_production_options(options)
    weights_trace = thread_trace === nothing ? nothing : thread_trace.weights
    lesser_trace = thread_trace === nothing ? nothing : thread_trace.lesser

    zeroΣ = zeros(ComplexF64, shape)
    η = effective_seed_parameters(problem).effective_scaled
    Gᴿ, κD, scaleD = _retarded_green_production(problem.grids.ε, h, zeroΣ; η = η, options)
    q = _seed_number_weights_production(
        Gᴿ,
        problem.grids,
        problem.physical.g_s;
        options,
        thread_trace = weights_trace,
    )

    kBT = _electronvolts(CODATA.kᴮₑᵥ * problem.physical.Tᴸ) / problem.scales.E₀_eV
    target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    number_at(μ) = _seed_number_from_weights(μ, problem.grids.ε, q, kBT) - target
    μlo = problem.grids.ε[1] - 50kBT
    μhi = problem.grids.ε[end] + 50kBT
    flo, fhi = number_at(μlo), number_at(μhi)
    flo * fhi ≤ 0 || throw(ArgumentError("seed chemical-potential root is not bracketed"))
    μ₀ = find_zero(number_at, (μlo, μhi), Bisection())

    Gˡ = zeros(ComplexF64, shape)
    _seed_lesser_production!(
        Gˡ,
        Gᴿ,
        problem.grids.ε,
        μ₀,
        kBT;
        options,
        thread_trace = lesser_trace,
    )
    # Bisection enforces the number at machine accuracy. A unilateral rescale
    # after this root would spoil exact Fermi bounds in nearly full tail blocks.
    return _green_state_production(Gᴿ, Gˡ, κD, scaleD, options), μ₀
end
