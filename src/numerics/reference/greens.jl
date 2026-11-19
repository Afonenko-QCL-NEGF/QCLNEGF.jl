_matrix_block(X, e, m) = Matrix(view(X, e, m, :, :))

"""
    retarded_green(ε, h, Σᴿ; η=0, return_scale=false)

Обратить полную плотную матрицу Дайсона для каждой пары `(E_e,k_m)`. Перед
`inv` применяется точное скалярное масштабирование по сумме в строке.
Оно сохраняет число обусловленности. Возвращаемая безразмерная матрица
`Gᴿ` имеет размер `(N_E,N_k,N_b,N_b)`.

Результат по умолчанию имеет вид `(Gᴿ, condition_number)`.
При `return_scale=true` третьим результатом служит безразмерный масштаб
уравновешивания ``s_{em}``.

Реализует [EQ-DYSON-001](@ref eq-dyson-retarded). Соглашения об индексах
и знаках описаны в [главе о функциях Грина и плотности](@ref theory-greens).
"""
function retarded_green(
    ε::AbstractVector{<:Real},
    h::AbstractArray{<:Number,3},
    Σᴿ::AbstractArray{<:Number,4};
    η::Real = 0.0,
    return_scale::Bool = false,
)
    NE, Nk, Nb, Nb2 = size(Σᴿ)
    Nb == Nb2 || throw(DimensionMismatch("Σᴿ blocks must be square"))
    size(h) == (Nk, Nb, Nb) || throw(DimensionMismatch("h has wrong shape"))
    length(ε) == NE || throw(DimensionMismatch("energy axis differs"))
    Gᴿ = zeros(ComplexF64, NE, Nk, Nb, Nb)
    κD = zeros(Float64, NE, Nk)
    scaleD = zeros(Float64, NE, Nk)
    I_b = Matrix{ComplexF64}(I, Nb, Nb)
    for e = 1:NE, m = 1:Nk
        D = (ε[e] + im * η) * I_b - Matrix(view(h, m, :, :)) - _matrix_block(Σᴿ, e, m)
        scale = maximum(sum(abs, D; dims = 2))
        scale > 0 || throw(SingularException(0))
        Gᴿ[e, m, :, :] .= inv(D / scale) / scale
        κD[e, m] = cond(D)
        scaleD[e, m] = scale
    end
    return return_scale ? (Gᴿ, κD, scaleD) : (Gᴿ, κD)
end

"""
Вычислить `A=i(Gᴿ-Gᴬ)` для каждого блока энергии и импульса.

См. [функции Грина и плотность](@ref theory-greens)
и [тождества Келдыша](@ref eq-keldysh).
"""
function spectral_function(Gᴿ::AbstractArray{<:Number,4})
    A = similar(ComplexF64.(Gᴿ))
    for e in axes(Gᴿ, 1), m in axes(Gᴿ, 2)
        GR = _matrix_block(Gᴿ, e, m)
        A[e, m, :, :] .= im .* (GR - GR')
    end
    return A
end

"""
    keldysh_green(Gᴿ, Σˡ)

Независимо вычислить полное плотное произведение Келдыша
`Gˡ=GᴿΣˡGᴬ` для каждой пары `(E_e,k_m)`.
Реализует [EQ-KELDYSH-001](@ref eq-keldysh).
См. [функции Грина и плотность](@ref theory-greens).
"""
function keldysh_green(Gᴿ::AbstractArray{<:Number,4}, Σˡ::AbstractArray{<:Number,4})
    size(Gᴿ) == size(Σˡ) || throw(DimensionMismatch("Gᴿ and Σˡ differ"))
    Gˡ = zeros(ComplexF64, size(Gᴿ))
    for e in axes(Gᴿ, 1), m in axes(Gᴿ, 2)
        GR = _matrix_block(Gᴿ, e, m)
        Gˡ[e, m, :, :] .= GR * _matrix_block(Σˡ, e, m) * GR'
    end
    return Gˡ
end

"""
Построить `Gᵍ=Gᴿ-Gᴬ+Gˡ` без отдельного хранения опережающей компоненты.

См. [функции Грина и плотность](@ref theory-greens)
и [тождества Келдыша](@ref eq-keldysh).
"""
function greater_green(Gᴿ::AbstractArray{<:Number,4}, Gˡ::AbstractArray{<:Number,4})
    size(Gᴿ) == size(Gˡ) || throw(DimensionMismatch("G arrays differ"))
    Gᵍ = zeros(ComplexF64, size(Gᴿ))
    for e in axes(Gᴿ, 1), m in axes(Gᴿ, 2)
        GR = _matrix_block(Gᴿ, e, m)
        Gᵍ[e, m, :, :] .= GR - GR' + _matrix_block(Gˡ, e, m)
    end
    return Gᵍ
end

"""
Вычислить безразмерную поверхностную концентрацию по меньшей функции
Грина с индексами `(E,k,a,b)`.

Реализует [EQ-NUMBER-001](@ref eq-number-normalization).
См. [функции Грина и плотность](@ref theory-greens).
"""
function number_functional(Gˡ::AbstractArray{<:Number,4}, grids::ModelGrids, g_s::Integer)
    NE, Nk, _, _ = size(Gˡ)
    value = 0.0 + 0.0im
    for e = 1:NE, m = 1:Nk
        value += grids.wᴱ[e] * grids.wᵏ[m] * tr(_matrix_block(Gˡ, e, m))
    end
    value *= -im * g_s / (2π)
    abs(imag(value)) ≤ _NUMBER_FUNCTIONAL_REALNESS_LIMIT * max(abs(real(value)), 1.0) ||
        throw(DomainError(value, "number functional is not real"))
    return real(value)
end

"""
    normalize_lesser(Gˡ, target, grids, g_s)

Односторонняя нормировка числа частиц для явно выбранного сравнительного
режима. Возвращает `(Gˡ_N, λ)` при `λ=N[Gˡ]/target`. Операция может нарушить
верхнюю границу Паули. В физической неподвижной точке выполняется `λ→1`.
По умолчанию SCBA использует парное преобразование `_normalize_keldysh_pair!`.

Реализует [EQ-NUMBER-001](@ref eq-number-normalization).
См. [функции Грина и плотность](@ref theory-greens).
"""
function normalize_lesser(
    Gˡ::AbstractArray{<:Number,4},
    target::Real,
    grids::ModelGrids,
    g_s::Integer,
)
    isfinite(target) && target > 0 ||
        throw(ArgumentError("target charge must be finite and positive"))
    λ = number_functional(Gˡ, grids, g_s) / target
    isfinite(λ) && λ > 0 || throw(DomainError(λ, "kinetic eigenvalue λ must be positive"))
    return ComplexF64.(Gˡ ./ λ), λ
end

function _green_state(Gᴿ, Gˡ, κD, scaleD)
    Gᵍ = greater_green(Gᴿ, Gˡ)
    A = spectral_function(Gᴿ)
    return GreenState(ComplexF64.(Gᴿ), ComplexF64.(Gˡ), Gᵍ, A, κD, scaleD)
end

"""
Коэффициенты преобразования Келдыша с заданным числом частиц и сохранением
ограничения Паули.

`N` и `H` — интегралы положительных матриц занятых и пустых состояний
`Gn=-iG<` и `Gp=iG>`. Избыток электронов переносится из `Gn` в `Gp`
в доле `(1-target/N)`. При дефиците доля `(target-N)/H` переносится из
`Gp` в `Gn`. Оба выпуклых преобразования сохраняют `Gn+Gp` и конус Паули
для допустимого исходного кандидата. Собственные значения не обрезаются.
При совпадении сырого числа частиц с целевым преобразование тождественно.
Сырое отношение `lambda=N/target` проверяется отдельно. Это ограничение
итерации. Оно не вводит распределение Ферми или физический резервуар.
"""
function _paired_number_coefficients(number::Real, holes::Real, target::Real)
    isfinite(target) && target > 0 ||
        throw(ArgumentError("target charge must be finite and positive"))
    isfinite(number) && number > 0 ||
        throw(DomainError(number, "raw occupied number must be finite and positive"))
    isfinite(holes) && holes >= 0 ||
        throw(DomainError(holes, "raw empty-state capacity must be finite and nonnegative"))
    target <= number + holes || throw(
        DomainError(
            (;
                target,
                raw_occupied = number,
                raw_unoccupied = holes,
                capacity = number+holes,
            ),
            "target charge exceeds the represented Pauli capacity (solver units); refine the domain/basis",
        ),
    )
    if number >= target
        return (;
            occupied_fraction = target/number,
            empty_fraction = 0.0,
            lambda = number/target,
        )
    end
    return (;
        occupied_fraction = 1.0,
        empty_fraction = (target-number)/holes,
        lambda = number/target,
    )
end

"""Применить парное ограничение числа частиц на месте к независимо вычисленным `G<` и `G>`."""
function _normalize_keldysh_pair!(Gˡ, Gᵍ, number::Real, holes::Real, target::Real)
    size(Gˡ) == size(Gᵍ) || throw(DimensionMismatch("Keldysh arrays differ"))
    weights = _paired_number_coefficients(number, holes, target)
    a, c = weights.occupied_fraction, weights.empty_fraction
    for index in eachindex(Gˡ, Gᵍ)
        lesser, greater = Gˡ[index], Gᵍ[index]
        # G< = i Gn and G> = -i Gp: transfers therefore have minus signs.
        Gˡ[index] = a*lesser - c*greater
        Gᵍ[index] = (1-c)*greater - (1-a)*lesser
    end
    return weights.lambda
end
