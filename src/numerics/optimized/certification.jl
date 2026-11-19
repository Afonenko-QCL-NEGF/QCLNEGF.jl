function _algorithm_options_from_solution(solution::NEGFSolution)
    contract = get(solution.observables, :restart_contract, solution.scba.restart_contract)
    contract === nothing && return AlgorithmOptions(
        solver_backend = :educational,
        energy_shift = :dense,
        hilbert = :direct,
        contraction = :literal,
        kernel_build = :direct,
        elastic_retarded = Symbol(get(solution.observables, :elastic_retarded, "hilbert")),
    )
    raw = get(contract, "algorithms", nothing)
    raw isa AbstractDict ||
        throw(ArgumentError("saved solver algorithm identity is unavailable"))
    defaults = AlgorithmOptions()
    return AlgorithmOptions(;
        (
            name=>(
                name === :elastic_retarded ? Symbol(get(raw, String(name), "hilbert")) :
                getfield(defaults, name) isa Symbol ? Symbol(raw[String(name)]) :
                raw[String(name)]
            ) for name in fieldnames(AlgorithmOptions)
        )...,
    )
end

"""Построить выбранное сырое физическое отображение для независимого итогового аудита."""
function _audit_candidate_for_solution(solution::NEGFSolution)
    algorithms = _algorithm_options_from_solution(solution)
    algorithms.solver_backend === :educational &&
        return _reference_audit_candidate(
            solution.problem, solution.scba; elastic_retarded = algorithms.elastic_retarded,
        )
    production = ProductionOptions(algorithms = algorithms)
    cache = build_production_cache(solution.problem; options = production)
    scattering, roundoff = _scattering_candidate_production(
        solution.problem,
        solution.scba.green,
        cache,
        production;
        return_roundoff = true,
    )
    embedding, plus, minus = _embedding_self_energy_selected(
        solution.problem,
        solution.scba.green,
        cache,
        production;
        scattering,
        h = project_hamiltonians(solution.problem, solution.Uᴴ),
    )
    return FixedPointAuditCandidate(scattering, embedding, plus, minus, roundoff)
end

"""
    certify_solution(solution; max_scba=2000, max_poisson=10, backend=:production)

Запустить ограниченную проверку с сохранённого состояния. Дополнительный
лимит итераций выделяется отдельно. Приближённое завершение и ранняя
остановка по стагнации отключены. Пороги не слабее исходных. При неудаче
возвращаются данные с результатом научной проверки. Физическая модель
и дискретизация сохраняются.
"""
function certify_solution(
    solution::NEGFSolution;
    max_scba::Integer = 2000,
    max_poisson::Integer = 10,
    backend::Symbol = :production,
    production_options::Union{Nothing,ProductionOptions} = nothing,
)
    max_scba > 0 && max_poisson > 0 ||
        throw(ArgumentError("certification budgets must be positive"))
    backend in (:production, :educational) ||
        throw(ArgumentError("unknown certification backend"))
    base = SolverTolerances()
    previous = solution.options.tolerances
    tolerances = SolverTolerances(;
        (
            name=>min(getfield(base, name), getfield(previous, name)) for
            name in fieldnames(SolverTolerances)
        )...,
    )
    old_policy = solution.options.convergence
    policy = ConvergencePolicy(;
        (
            name=>(
                name === :mode ? :strict_fail_fast :
                name === :diagnostic_quality ?
                DiagnosticQualityPolicy(enabled = false) :
                name === :stagnation_window ? 0 :
                name === :stagnation_relative_improvement ? 0.0 :
                getfield(old_policy, name)
            ) for name in fieldnames(ConvergencePolicy)
        )...,
    )
    old = solution.options
    options = SolverOptions(;
        (
            name=>(
                name === :max_scba ? Int(max_scba) :
                name === :max_poisson ? Int(max_poisson) :
                name === :tolerances ? tolerances :
                name === :convergence ? policy : getfield(old, name)
            ) for name in fieldnames(SolverOptions)
        )...,
    )
    if backend === :educational
        _algorithm_options_from_solution(solution).solver_backend === :educational || throw(
            ArgumentError("educational certification cannot change the saved physical map"),
        )
        result = solve(
            solution.problem;
            options,
            initial_Uᴴ = solution.Uᴴ,
            initial_scba = solution.scba,
            elastic_retarded = _algorithm_options_from_solution(solution).elastic_retarded,
        )
    else
        algorithms = _algorithm_options_from_solution(solution)
        production =
            production_options === nothing ? ProductionOptions(algorithms = algorithms) :
            production_options
        _restart_contract_value(production.algorithms) ==
        _restart_contract_value(algorithms) || throw(
            ArgumentError(
                "certification algorithm identity differs from the published state",
            ),
        )
        result = _solve_production(
            solution.problem;
            options,
            production_options = production,
            initial_Uᴴ = solution.Uᴴ,
            initial_scba = solution.scba,
        )
    end
    result.observables[:certification] = Dict{String,Any}(
        "explicit"=>true,
        "max_scba"=>Int(max_scba),
        "max_poisson"=>Int(max_poisson),
        "source_status"=>String(solution.status),
        "passed"=>result.converged,
        "policy"=>"strict_fail_fast",
        "approximate_allowed"=>false,
    )
    return result
end
