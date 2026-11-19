function _fermi(ε::Real, μ::Real, kBT::Real)
    x = (ε - μ) / kBT
    x > _FERMI_EXPONENT_SATURATION && return 0.0
    x < -_FERMI_EXPONENT_SATURATION && return 1.0
    return inv(exp(x) + 1)
end

"""
Запрошенное и фактическое начальное уширение в единицах решателя и в эВ.

Запрошенное значение применяется точно. Шаг сетки сообщается
для проверки разрешения. При уточнении сетки начальное уширение
остаётся заданным. В ходе итераций оно не применяется.
"""
function effective_seed_parameters(problem::NEGFProblem)
    return _effective_seed_parameters(
        problem.numerical,
        problem.scales,
        problem.grids.wᴱ[2],
    )
end

function effective_seed_parameters(
    numerical::NumericalParameters,
    scales::ScaleSystem = ScaleSystem(),
)
    step = scaled_value(numerical.E_max-numerical.E_min, scales, :energy)/(numerical.N_E-1)
    grid_floor = numerical.N_E == 2 ? step/2 : step
    return _effective_seed_parameters(numerical, scales, grid_floor)
end

function _effective_seed_parameters(
    numerical::NumericalParameters,
    scales::ScaleSystem,
    grid_floor::Real,
)
    requested = scaled_value(numerical.η_seed, scales, :energy)
    effective = requested
    scale = scales.E₀_eV
    return (;
        requested_scaled = requested,
        effective_scaled = effective,
        grid_floor_scaled = grid_floor,
        requested_eV = requested*scale,
        effective_eV = effective*scale,
        grid_floor_eV = grid_floor*scale,
        clamped = effective > requested,
        rule = :explicit_requested,
        working_eta_eV = 0.0,
    )
end

function _seed_green(problem::NEGFProblem, h::AbstractArray{<:Number,3})
    n = problem.numerical
    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    zeroΣ = zeros(ComplexF64, shape)
    η = effective_seed_parameters(problem).effective_scaled
    Gᴿ, κD, scaleD = retarded_green(problem.grids.ε, h, zeroΣ; η = η, return_scale = true)
    A = spectral_function(Gᴿ)
    kBT = _electronvolts(CODATA.kᴮₑᵥ * problem.physical.Tᴸ) / problem.scales.E₀_eV
    target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ

    function number_at(μ)
        Gˡ = zeros(ComplexF64, shape)
        for e = 1:n.N_E
            Gˡ[e, :, :, :] .= im * _fermi(problem.grids.ε[e], μ, kBT) .* A[e, :, :, :]
        end
        return number_functional(Gˡ, problem.grids, problem.physical.g_s) - target
    end
    μlo = problem.grids.ε[1] - 50kBT
    μhi = problem.grids.ε[end] + 50kBT
    flo, fhi = number_at(μlo), number_at(μhi)
    flo * fhi ≤ 0 || throw(ArgumentError("seed chemical-potential root is not bracketed"))
    μ₀ = find_zero(number_at, (μlo, μhi), Bisection())
    Gˡ = zeros(ComplexF64, shape)
    for e = 1:n.N_E
        Gˡ[e, :, :, :] .= im * _fermi(problem.grids.ε[e], μ₀, kBT) .* A[e, :, :, :]
    end
    return _green_state(Gᴿ, Gˡ, κD, scaleD), μ₀
end

function _copy_family(f::SelfEnergyFamily)
    return SelfEnergyFamily(copy(f.Σᴿ), copy(f.Σˡ), copy(f.Σᵍ))
end

function _copy_mechanisms(d::Dict{Symbol,SelfEnergyFamily})
    return Dict(name => _copy_family(family) for (name, family) in d)
end

function _total_family(
    scattering::Dict{Symbol,SelfEnergyFamily},
    embedding::SelfEnergyFamily,
    shape::NTuple{4,Int},
)
    sc = _sum_selfenergies(scattering, shape)
    return SelfEnergyFamily(
        sc.Σᴿ + embedding.Σᴿ,
        sc.Σˡ + embedding.Σˡ,
        sc.Σᵍ + embedding.Σᵍ,
    )
end

function _current_green(problem::NEGFProblem, h, scattering, embedding)
    shape = (
        problem.numerical.N_E,
        problem.numerical.N_k,
        problem.numerical.N_b,
        problem.numerical.N_b,
    )
    total = _total_family(scattering, embedding, shape)
    Gᴿ, κD, scaleD =
        retarded_green(problem.grids.ε, h, total.Σᴿ; η = 0.0, return_scale = true)
    raw = keldysh_green(Gᴿ, total.Σˡ)
    greater = keldysh_green(Gᴿ, total.Σᵍ)
    target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    number = number_functional(raw, problem.grids, problem.physical.g_s)
    holes = -number_functional(greater, problem.grids, problem.physical.g_s)
    λ = _normalize_keldysh_pair!(raw, greater, number, holes, target)
    return GreenState(Gᴿ, raw, greater, spectral_function(Gᴿ), κD, scaleD), λ, total
end

function _combined_candidate(
    scattering::Dict{Symbol,SelfEnergyFamily},
    embedding::SelfEnergyFamily,
    shape::NTuple{4,Int},
)
    return _total_family(scattering, embedding, shape)
end

"""
    solve_scba(problem, Uᴴ; options=baseline_options(), initial=nothing)

Решить бесконтактную задачу Дайсона — Келдыша — SCBA при заданном числе
частиц и фиксированной безразмерной энергии Хартри. Все кандидатные
собственные энергии строятся отображением Якоби из одной общей
`G_ν`. Начальное уширение полностью удаляется (`η_work=0`).

Реализует [внутренний цикл SCBA](@ref theory-scba).
См. [функции Грина и плотность](@ref theory-greens)
и [микроскопические ядра](@ref theory-kernels).
"""
function solve_scba(
    problem::NEGFProblem,
    Uᴴ::AbstractVector{<:Real};
    options::SolverOptions = baseline_options(),
    initial::Union{Nothing,SCBAResult} = nothing,
    history_observer::Union{Nothing,Function} = nothing,
    elastic_retarded::Symbol = :hilbert,
)
    _check_options(options)
    elastic_retarded in (:hilbert, :direct_green) ||
        throw(ArgumentError("elastic_retarded must be :hilbert or :direct_green"))
    reference_algorithms = AlgorithmOptions(;
        solver_backend = :educational, energy_shift = :dense,
        hilbert = :direct, contraction = :literal, kernel_build = :direct,
        elastic_retarded,
    )
    # Preserve the historical default result; explicitly tag the new variant.
    restart_contract = elastic_retarded === :hilbert ? nothing :
                       _solver_restart_contract(options, reference_algorithms)
    if initial !== nothing
        saved_mode = initial.restart_contract === nothing ? "hilbert" :
                     get(get(initial.restart_contract, "algorithms", Dict()),
                         "elastic_retarded", "hilbert")
        saved_mode == String(elastic_retarded) || throw(ArgumentError(
            "reference warm state elastic retarded closure differs; start a new run",
        ))
    end
    h = project_hamiltonians(problem, Uᴴ)
    n = problem.numerical
    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    history = SCBAIteration[]

    if initial === nothing
        green, _ = _seed_green(problem, h)
        scattering = _scattering_candidate(problem, green; elastic_retarded)
        embedding, plus, minus = embedding_self_energy(problem, green)
    else
        Set(keys(initial.scattering)) == Set(problem.kernels.enabled) ||
            throw(ArgumentError("initial SCBA mechanisms do not match problem"))
        all(
            size(getfield(family, component)) == shape for
            family in values(initial.scattering) for component in (:Σᴿ, :Σˡ, :Σᵍ)
        ) || throw(DimensionMismatch("initial scattering arrays have wrong shape"))
        all(
            size(getfield(family, component)) == shape for
            family in (initial.embedding, initial.embedding_plus, initial.embedding_minus)
            for component in (:Σᴿ, :Σˡ, :Σᵍ)
        ) || throw(DimensionMismatch("initial embedding arrays have wrong shape"))
        all(
            size(getfield(initial.green, component)) == shape for
            component in (:Gᴿ, :Gˡ, :Gᵍ, :A)
        ) || throw(DimensionMismatch("initial Green arrays have wrong shape"))
        size(initial.green.condition_number) == (n.N_E, n.N_k) &&
        size(initial.green.dyson_scale) == (n.N_E, n.N_k) ||
            throw(DimensionMismatch("initial Green diagnostics have wrong shape"))
        green = initial.green
        scattering = _copy_mechanisms(initial.scattering)
        embedding = _copy_family(initial.embedding)
        plus = _copy_family(initial.embedding_plus)
        minus = _copy_family(initial.embedding_minus)
    end
    Jprevious = nothing
    populations_previous = nothing
    convergence_streak = 0
    diagnostic_quality_streak = 0

    for ν = 1:options.max_scba
        # Do not convert a Dyson/Keldysh failure into an SCBAResult: `green`
        # still describes the warm state and is not aligned with this `h`.
        # The caller must retain the last atomically checkpointed state instead.
        green_current, λ, total = _current_green(problem, h, scattering, embedding)
        green = green_current

        local candidate_sc, candidate_embedding, candidate_plus, candidate_minus, rround
        try
            candidate_sc, rround =
                _scattering_candidate(problem, green; return_roundoff = true, elastic_retarded)
            candidate_embedding, candidate_plus, candidate_minus =
                embedding_self_energy(problem, green)
        catch error
            error isa DomainError || rethrow()
            @warn "SCBA stopped while building a candidate" iteration=ν exception=(
                error,
                catch_backtrace(),
            )
            _record_candidate_failure!(history, problem, h, green, total, plus, ν, λ)
            return SCBAResult(
                green,
                scattering,
                embedding,
                plus,
                minus,
                history,
                false,
                :invalid_candidate,
                :invalid,
                restart_contract,
            )
        end
        candidate_total = _combined_candidate(candidate_sc, candidate_embedding, shape)

        rD = _dyson_residual(problem.grids.ε, h, total.Σᴿ, green.Gᴿ)
        rA = _spectral_residual(green, total)
        rK = _keldysh_residual(green, candidate_total)
        rΣ = _selfenergy_residual(embedding, candidate_embedding)
        rΣ = max(
            rΣ,
            _selfenergy_residual(plus, candidate_plus),
            _selfenergy_residual(minus, candidate_minus),
        )
        for mechanism in keys(candidate_sc)
            rΣ = max(
                rΣ,
                _selfenergy_residual(scattering[mechanism], candidate_sc[mechanism]),
            )
        end
        rλ = abs(λ - 1)
        rPSD = _green_positivity(green, total)
        rcaus = _causality_residual(total)
        Jbar = _boundary_flux_bar(problem, plus, green)
        J = Jbar * Float64(ustrip(u"A/m^2", problem.scales.J₀))
        target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
        populations = real.(diag(_sheet_density_matrix_bar(problem, green.Gˡ))) ./ target
        rJchange =
            Jprevious === nothing ? Inf :
            abs(J - Jprevious) / max(abs(J), abs(Jprevious), 1e-12)
        rpopulation =
            populations_previous === nothing ? Inf :
            norm(populations - populations_previous) / (norm(populations) + 1e-14)
        push!(
            history,
            SCBAIteration(
                ν,
                rD,
                rA,
                rK,
                rΣ,
                rλ,
                λ,
                rPSD,
                rcaus,
                rround,
                rJchange,
                rpopulation,
                J,
                λ*target,
                target,
                number_functional(green.Gˡ, problem.grids, problem.physical.g_s),
                isempty(history) ? NaN : λ-last(history).λ,
                nothing,
                nothing,
            ),
        )

        history_observer===nothing || history_observer(last(history))
        if !all(isfinite, (rD, rA, rK, rΣ, rλ, rPSD, rcaus, rround, J, λ))
            return SCBAResult(
                green,
                scattering,
                embedding,
                plus,
                minus,
                history,
                false,
                :nonfinite_metrics,
                :invalid,
                restart_contract,
            )
        end
        assessment = scba_convergence_assessment(
            last(history),
            options.tolerances,
            options.convergence,
        )
        convergence_streak = assessment.passed ? convergence_streak + 1 : 0
        diagnostic = scba_diagnostic_quality_assessment(
            last(history),
            options.tolerances,
            options.convergence,
        )
        diagnostic_state = _scba_diagnostic_quality_state(
            diagnostic_quality_streak,
            diagnostic,
            options.convergence,
        )
        diagnostic_quality_streak = diagnostic_state.streak
        quality = diagnostic_state.quality
        decision = scba_iteration_decision(history, options, convergence_streak)
        if decision.terminal
            return SCBAResult(
                green,
                scattering,
                embedding,
                plus,
                minus,
                history,
                decision.converged,
                decision.status,
                decision.quality,
                restart_contract,
            )
        end

        scattering = _mix_mechanisms(scattering, candidate_sc, options.α_Σ)
        plus = _mix(plus, candidate_plus, options.α_Σ)
        minus = _mix(minus, candidate_minus, options.α_Σ)
        embedding =
            SelfEnergyFamily(plus.Σᴿ + minus.Σᴿ, plus.Σˡ + minus.Σˡ, plus.Σᵍ + minus.Σᵍ)
        Jprevious = J
        populations_previous = populations
    end
    error("unreachable SCBA loop exit")
end

"""Сохранить попытку итерации, если для неё не определён физический кандидат."""
function _record_candidate_failure!(
    history,
    problem,
    h,
    green,
    total,
    plus,
    iteration,
    lambda,
)
    target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    current =
        _boundary_flux_bar(problem, plus, green)*Float64(
            ustrip(u"A/m^2", problem.scales.J₀),
        )
    push!(
        history,
        SCBAIteration(
            iteration,
            _dyson_residual(problem.grids.ε, h, total.Σᴿ, green.Gᴿ),
            _spectral_residual(green, total),
            Inf,
            Inf,
            abs(lambda-1),
            lambda,
            _green_positivity(green, total),
            _causality_residual(total),
            Inf,
            Inf,
            Inf,
            current,
            lambda*target,
            target,
            number_functional(green.Gˡ, problem.grids, problem.physical.g_s),
            isempty(history) ? NaN : lambda-last(history).λ,
            nothing,
            nothing,
        ),
    )
    return history
end
