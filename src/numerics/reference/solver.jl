"""
    build_problem(; physical=reference_parameters(), numerical=baseline_numerics(),
                    scattering=default_scattering(), scales=ScaleSystem())

Build all immutable data: grids, exactly normalized profiles, fixed localized
basis, complete kernels, and four dense non-periodic shift matrices.  This is
the expensive static stage of the direct reference algorithm.

See [the implementation plan](@ref implementation-plan),
[the array contract](@ref array-contracts), and
[static verification](@ref theory-validation).
"""
function build_problem(;
    physical::PhysicalParameters = reference_parameters(),
    numerical::NumericalParameters = baseline_numerics(),
    scattering::ScatteringOptions = default_scattering(),
    scales::ScaleSystem = ScaleSystem(),
    localization::Symbol = :pzp,
    validate_static::Bool = true,
    physical_models::PhysicalModelOptions = PhysicalModelOptions(),
)
    physical = resolve_physical_models(physical, physical_models)
    grids = build_grids(physical, numerical, scales)
    profiles = build_profiles(physical, numerical, scales, grids)
    basis = build_basis(physical, numerical, scales, grids, profiles; localization)
    kernels = build_kernels(physical, numerical, scattering, scales, grids, profiles, basis)
    sp = _scaled_physics(physical, scales)
    W₊ᴱᵖ = build_shift_matrix(grids.ε, +sp.Eᵖ)
    W₋ᴱᵖ = build_shift_matrix(grids.ε, -sp.Eᵖ)
    W₊ᴸᴼ = build_shift_matrix(grids.ε, +sp.ħωᴸᴼ)
    W₋ᴸᴼ = build_shift_matrix(grids.ε, -sp.ħωᴸᴼ)
    problem = NEGFProblem(
        physical,
        numerical,
        scattering,
        scales,
        grids,
        profiles,
        basis,
        kernels,
        W₊ᴱᵖ,
        W₋ᴱᵖ,
        W₊ᴸᴼ,
        W₋ᴸᴼ,
    )
    augmented_kernels = add_sppa_kernel(
        problem.physical,
        problem.scales,
        problem.grids,
        problem.basis,
        problem.kernels,
        physical_models.electron_electron,
    )
    problem = NEGFProblem(
        problem.physical,
        problem.numerical,
        problem.scattering,
        problem.scales,
        problem.grids,
        problem.profiles,
        problem.basis,
        augmented_kernels,
        problem.W₊ᴱᵖ,
        problem.W₋ᴱᵖ,
        problem.W₊ᴸᴼ,
        problem.W₋ᴸᴼ,
        problem.energy_shift_discretization,
        physical_models,
    )

    if validate_static
        report = validate_problem(problem)
        report.passed || throw(
            ArgumentError("static model validation failed: " * join(report.messages, "; ")),
        )
    end
    return problem
end

function _outer_current_residual(problem, scba)
    plus = _boundary_flux_bar(problem, scba.embedding_plus, scba.green)
    minus = _boundary_flux_bar(problem, scba.embedding_minus, scba.green)
    residual = abs(plus + minus) / max(abs(plus), abs(minus), 1e-12)
    J = plus * Float64(ustrip(u"A/m^2", problem.scales.J₀))
    return residual, J
end

function _collect_observables(
    problem::NEGFProblem,
    scba::SCBAResult,
    Uᴴ::AbstractVector{<:Real},
)
    green = scba.green
    populations = state_populations(problem, green.Gˡ)
    current = CODATA.e * boundary_flux(problem, scba.embedding_plus, green)
    return Dict{Symbol,Any}(
        :sheet_density_matrix => sheet_density_matrix(problem, green.Gˡ),
        :state_populations => populations,
        :effective_levels => effective_levels(problem, Uᴴ),
        :electron_flow_current => uconvert(u"A/m^2", current),
        :energy_resolved_current =>
            energy_resolved_current(problem, scba.embedding_plus, green),
        :bond_current => bond_current(problem, green.Gˡ),
        :power_balance =>
            power_balance(problem, scba.scattering, scba.embedding_plus, green),
        :spectral_maps => spectral_maps(problem, green),
    )
end

"""
    solve(problem; options=baseline_options())

Run the nested fixed-density SCBA and periodic Hartree–Poisson loops.  The
returned `NEGFSolution` always carries an explicit status; reaching an
iteration limit is never reported as convergence. The last complete SCBA state
is published with its own accepted `Uᴴ`; a fresh candidate audit certifies the
same operator/state without an implicit unbudgeted extra solve.

Implements the [closed Poisson–SCBA algorithm](@ref theory-outer-loop).
Its inner map and final acceptance gates are detailed in
[the SCBA chapter](@ref theory-scba) and
[Verification](@ref theory-validation).
"""
function solve(
    problem::NEGFProblem;
    options::SolverOptions = baseline_options(),
    initial_Uᴴ::Union{Nothing,AbstractVector{<:Real}} = nothing,
    initial_scba::Union{Nothing,SCBAResult} = nothing,
    history_observer::Union{Nothing,Function} = nothing,
)
    _check_options(options)
    Nz = problem.numerical.N_z
    Uᴴ = initial_Uᴴ === nothing ? zeros(Float64, Nz) : Float64.(initial_Uᴴ)
    length(Uᴴ) == Nz && all(isfinite, Uᴴ) ||
        throw(ArgumentError("invalid initial Hartree field"))
    published_U = copy(Uᴴ)
    nprevious = nothing
    previous_scba = initial_scba
    Jprevious = nothing
    populations_previous = nothing
    outer_history = OuterIteration[]
    warnings = Dict{String,Any}[]
    outer_converged = false
    status = :max_poisson_iterations
    convergence_streak = 0

    for μ = 1:options.max_poisson
        inner_policy = inner_working_policy(options, outer_history)
        scba = solve_scba(
            problem,
            Uᴴ;
            options = inner_policy.options,
            initial = previous_scba,
            history_observer = history_observer===nothing ? nothing :
                               (row->history_observer(:scba, μ, row, problem)),
        )
        published_U = copy(Uᴴ)
        if options.convergence.mode === :research_continue && !scba.converged
            push!(warnings, _research_transition_warning(scba, options, μ))
        end
        if !scba_accepted(scba)
            previous_scba = scba
            status = Symbol(:scba_, scba.status)
            break
        end
        if scba.status === :approximate
            push!(warnings, _approximate_warning(scba, inner_policy.options, μ))
        end
        n̄ = _electron_density_bar(problem, scba.green.Gˡ)
        candidate, ζ, rP, rneutral = solve_periodic_poisson(problem, n̄)
        Unew = (1 - options.α_P) .* Uᴴ .+ options.α_P .* candidate
        rU = maximum(abs.(candidate .- Uᴴ))
        rn =
            nprevious === nothing ? Inf :
            sqrt(sum(problem.grids.wˣ .* (n̄ .- nprevious) .^ 2)) / (
                _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ /
                sqrt(_scaled_physics(problem.physical, problem.scales).Lp)
            )
        rJ, J = _outer_current_residual(problem, scba)
        populations =
            real.(diag(_sheet_density_matrix_bar(problem, scba.green.Gˡ))) ./
            _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
        rJchange =
            Jprevious === nothing ? Inf :
            abs(J - Jprevious) / max(abs(J), abs(Jprevious), 1e-12)
        rpopulation =
            populations_previous === nothing ? Inf :
            norm(populations - populations_previous) / (norm(populations) + 1e-14)
        push!(
            outer_history,
            OuterIteration(
                μ,
                rP,
                rU,
                rn,
                rneutral,
                rJ,
                rJchange,
                rpopulation,
                ζ,
                J,
                copy(populations),
                copy(n̄),
            ),
        )
        history_observer===nothing ||
            history_observer(:outer, μ, last(outer_history), problem)
        Uᴴ = Unew
        nprevious = n̄
        previous_scba = scba
        Jprevious = J
        populations_previous = populations
        assessment = poisson_convergence_assessment(
            last(outer_history),
            options.tolerances,
            options.convergence,
        )
        convergence_streak =
            assessment.passed && scba.converged ? convergence_streak + 1 : 0
        if convergence_streak >= options.convergence.required_consecutive_poisson_passes
            outer_converged = true
            status = :converged
            break
        end
    end

    # Publish the last complete fixed-Hartree snapshot. The unexecuted next
    # Hartree update is not a solution and does not cause an implicit extra solve.
    final_scba = previous_scba
    Uᴴ = published_U
    nfinal = _electron_density_bar(problem, final_scba.green.Gˡ)
    observables = _collect_observables(problem, final_scba, Uᴴ)
    observables[:warnings] = warnings
    observables[:acceptance_measurement_context] = Dict{String,Any}(
        "backend"=>"educational",
        "hilbert_roundoff_measured"=>true,
        "fresh_candidate_required"=>true,
    )
    dummy = ConvergenceReport(false, Dict{Symbol,Float64}(), String[])
    provisional = NEGFSolution(
        problem,
        options,
        Uᴴ,
        nfinal,
        final_scba,
        outer_history,
        observables,
        dummy,
        false,
        status,
    )
    quality_report = _final_quality_report(provisional, outer_converged)
    _publish_charge_diagnostics!(observables, quality_report.report)
    _publish_model_diagnostics!(observables, problem, final_scba, Uᴴ)
    final_warning = _final_validation_warning(provisional, quality_report)
    final_warning === nothing || push!(warnings, final_warning)
    report = quality_report.report
    converged = quality_report.converged
    final_status =
        converged ? :converged :
        quality_report.approximate ? :approximate :
        (
            outer_converged && !final_scba.converged ? :final_scba_failed :
            outer_converged && !report.passed ? :validation_failed : status
        )
    result = NEGFSolution(
        problem,
        options,
        Uᴴ,
        nfinal,
        final_scba,
        outer_history,
        observables,
        report,
        converged,
        final_status,
    )
    observables[:quality] = solution_quality(result)
    observables[:termination_reason] = String(final_status)
    return result
end

function solve(;
    physical::PhysicalParameters = reference_parameters(),
    numerical::NumericalParameters = baseline_numerics(),
    scattering::ScatteringOptions = default_scattering(),
    scales::ScaleSystem = ScaleSystem(),
    options::SolverOptions = baseline_options(),
)
    problem = build_problem(; physical, numerical, scattering, scales)
    return solve(problem; options)
end
