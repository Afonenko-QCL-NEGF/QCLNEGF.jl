"""Evaluate the reference raw map; no optimized kernel is used as its oracle."""
function _reference_audit_candidate(problem::NEGFProblem, scba::SCBAResult)
    scattering, roundoff =
        _scattering_candidate(problem, scba.green; return_roundoff = true)
    embedding, plus, minus = embedding_self_energy(problem, scba.green)
    return FixedPointAuditCandidate(scattering, embedding, plus, minus, roundoff)
end

function charge_constraint_diagnostics(
    raw_charge::Real,
    target_charge::Real,
    normalized_charge::Real;
    previous_lambda::Real = NaN,
)
    all(isfinite, (raw_charge, target_charge, normalized_charge)) ||
        throw(ArgumentError("charge diagnostics require finite particle numbers"))
    target_charge > 0 && raw_charge > 0 ||
        throw(ArgumentError("raw and target particle numbers must be positive"))
    lambda = raw_charge/target_charge
    return ChargeConstraintDiagnostics(
        Float64(raw_charge),
        Float64(target_charge),
        Float64(normalized_charge),
        Float64(lambda),
        isfinite(previous_lambda) ? Float64(lambda-previous_lambda) : NaN,
        Float64(abs(lambda-1)),
        Float64(abs(normalized_charge-target_charge)/target_charge),
    )
end

"""Fresh coupled kinetic residuals, evaluated before any normalization of a candidate.

The selected physical map supplies candidate self-energies. Residual matrix
algebra is evaluated independently of the production residual pass. In
particular r_K compares published G< to GR Σ_candidate< GA, without dividing
the latter by a fitted particle-number eigenvalue.
"""
function fresh_fixed_point_metrics(
    problem::NEGFProblem,
    scba::SCBAResult,
    U::AbstractVector{<:Real},
    candidate::FixedPointAuditCandidate,
)
    shape = size(scba.green.Gᴿ)
    total = _total_family(scba.scattering, scba.embedding, shape)
    candidate_total = _total_family(candidate.scattering, candidate.embedding, shape)
    h = project_hamiltonians(problem, U)
    raw = keldysh_green(scba.green.Gᴿ, total.Σˡ)
    raw_number = number_functional(raw, problem.grids, problem.physical.g_s)
    normalized = number_functional(scba.green.Gˡ, problem.grids, problem.physical.g_s)
    target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    previous_lambda = length(scba.history) > 1 ? scba.history[end-1].λ : NaN
    charge = charge_constraint_diagnostics(raw_number, target, normalized; previous_lambda)
    candidate_raw_number = number_functional(
        keldysh_green(scba.green.Gᴿ, candidate_total.Σˡ),
        problem.grids,
        problem.physical.g_s,
    )
    candidate_lambda = candidate_raw_number/target
    r_sigma = maximum((
        _selfenergy_residual(scba.embedding, candidate.embedding),
        _selfenergy_residual(scba.embedding_plus, candidate.plus),
        _selfenergy_residual(scba.embedding_minus, candidate.minus),
    ))
    for name in problem.kernels.enabled
        r_sigma = max(
            r_sigma,
            _selfenergy_residual(scba.scattering[name], candidate.scattering[name]),
        )
    end
    return Dict{Symbol,Float64}(
        :r_D=>_dyson_residual(problem.grids.ε, h, total.Σᴿ, scba.green.Gᴿ),
        :r_D_candidate=>_dyson_residual(
            problem.grids.ε,
            h,
            candidate_total.Σᴿ,
            scba.green.Gᴿ,
        ),
        :r_A_candidate=>_spectral_residual(scba.green, candidate_total),
        :r_A=>_spectral_residual(scba.green, total),
        :r_K=>_keldysh_residual(scba.green, candidate_total),
        :r_Σ=>r_sigma,
        :r_λ=>charge.constraint_residual,
        :r_roundoff=>candidate.roundoff,
        :raw_charge=>charge.raw_charge,
        :target_charge=>charge.target_charge,
        :normalized_charge=>charge.normalized_charge,
        :lambda=>charge.lambda,
        :lambda_change=>charge.lambda_change,
        :r_charge_normalized=>charge.normalization_residual,
        :candidate_raw_charge=>candidate_raw_number,
        :candidate_lambda=>candidate_lambda,
        :r_lambda_candidate=>abs(candidate_lambda-1),
    )
end

"""Store exact final constraint values, independent of the neutralized Poisson source."""
function _publish_charge_diagnostics!(observables, report)
    observables[:charge_constraint] = Dict(
        String(name)=>get(report.metrics, name, NaN) for name in (
            :raw_charge,
            :target_charge,
            :normalized_charge,
            :lambda,
            :lambda_change,
            :r_λ,
            :r_charge_normalized,
            :candidate_raw_charge,
            :candidate_lambda,
            :r_lambda_candidate,
        )
    )
    return observables
end

const _NONLINEAR_ACCEPTANCE_METRICS = Set((
    :r_K,
    :r_Σ,
    :r_λ,
    :r_lambda_candidate,
    :r_D_candidate,
    :r_A_candidate,
    :r_Jchange_scba,
    :r_population_scba,
    :r_U,
    :r_n,
    :r_Jchange,
    :r_population,
))

"""Separate fixed-point progress, physical gates and certified acceptance.

Thresholds come from the same evaluator used by `validate_solution`, never an
API/UI replica. Missing measurements are explicit and cannot become zero or
pass. All numeric thresholds are project acceptance policy; these numbers are
not attributed to universal NEGF identities or borrowed literature tolerances.
"""
function solution_scientific_assessment(solution::NEGFSolution)
    limits =
        _solution_validation_limits(solution.options.tolerances, solution.report.metrics)
    records = Dict{String,Any}()
    nonlinear_pass, physical_pass = true, true
    for (name, threshold) in limits
        measured = haskey(solution.report.metrics, name)
        value = get(solution.report.metrics, name, NaN)
        metadata = _acceptance_metric_metadata(name, solution)
        applicable=metadata["applicability"]["applies"]
        scope=metadata["applicability"]["scope"]
        missing_prerequisite=(scope=="outer_history" && isempty(solution.outer_history)) ||
                             (scope=="inner_history" && isempty(solution.scba.history))
        measurement_context=get(
            solution.observables,
            :acceptance_measurement_context,
            Dict(),
        )
        unmeasured_transform=name===:r_roundoff &&
                             get(
            measurement_context,
            "hilbert_roundoff_measured",
            nothing,
        )===false
        status =
            !applicable ? "not_applicable" :
            (!measured || missing_prerequisite || unmeasured_transform) ? "not_measured" :
            !isfinite(value) ? "error" : value <= threshold ? "pass" : "fail"
        category =
            name in _NONLINEAR_ACCEPTANCE_METRICS ? "nonlinear_fixed_point" :
            "physical_and_algebraic"
        if name in _NONLINEAR_ACCEPTANCE_METRICS
            nonlinear_pass &= status in ("pass", "not_applicable")
        else
            physical_pass &= status in ("pass", "not_applicable")
        end
        records[String(name)] = merge(
            metadata,
            Dict{String,Any}(
                "value"=>status in ("pass", "fail") ? value : nothing,
                "recorded_value"=>isfinite(value) ? value : nothing,
                "status"=>status,
                "threshold"=>threshold,
                "comparison"=>"less_or_equal",
                "iteration_comparison"=>"inner/outer iteration confirmation uses strictly_less; final validation uses less_or_equal",
                "category"=>category,
                "reason"=>!applicable ? metadata["applicability"]["reason"] :
                          missing_prerequisite ?
                          "Required iteration history is unavailable." :
                          unmeasured_transform ?
                          "The runtime did not perform the independent Hilbert roundoff comparison." :
                          !measured ? "metric was not measured" :
                          !isfinite(value) ?
                          "non-finite measurement or unavailable prerequisite" :
                          nothing,
            ),
        )
    end
    streak = _trailing_pass_count(
        solution.outer_history,
        row -> poisson_convergence_assessment(
            row,
            solution.options.tolerances,
            solution.options.convergence,
        ).passed,
    )
    iterative =
        nonlinear_pass &&
        solution.scba.converged &&
        streak >= solution.options.convergence.required_consecutive_poisson_passes
    return Dict{String,Any}(
        "registry_version"=>_ACCEPTANCE_METADATA_VERSION,
        "registry_source"=>"src/numerics/reference/acceptance_metadata.jl",
        "model_context"=>Dict(
            "boundary_condition"=>"stationary_field_periodic_embedding",
            "basis_contract"=>"orthonormal finite projected basis",
            "energy_shift_discretization"=>String(
                solution.problem.energy_shift_discretization,
            ),
            "enabled_scattering"=>String.(solution.problem.kernels.enabled),
            "physical_models"=>_restart_contract_value(solution.problem.models),
        ),
        "quadrature"=>Dict(
            "energy"=>Dict("weights"=>copy(solution.problem.grids.wᴱ), "units"=>"E0"),
            "radial"=>Dict(
                "weights"=>copy(solution.problem.grids.wᵏ),
                "units"=>"L0^-2; includes radial 1/(2*pi) measure",
            ),
            "spatial"=>Dict("weights"=>copy(solution.problem.grids.wˣ), "units"=>"L0"),
            "spin_degeneracy"=>solution.problem.physical.g_s,
        ),
        "scales"=>Dict(
            "E0_eV"=>solution.problem.scales.E₀_eV,
            "L0_m"=>solution.problem.scales.L₀_m,
            "lambda_P"=>solution.problem.scales.λ_P,
            "Nb"=>solution.problem.numerical.N_b,
            "energy_tail_window_fraction"=>solution.options.energy_tail_window_fraction,
            "momentum_tail_window_fraction"=>solution.options.momentum_tail_window_fraction,
        ),
        "iterative_converged"=>iterative,
        "fixed_hartree_converged"=>solution.scba.converged,
        "physical_gates_passed"=>physical_pass,
        "stationary_candidate_accepted"=>solution.converged && iterative && physical_pass,
        "discretization_verified"=>"not_measured",
        "scientific_accepted"=>false,
        "scientific_acceptance_reason"=>"campaign-level independent discretization evidence is required",
        "outer_confirmation_count"=>streak,
        "metrics"=>records,
        "evaluator"=>"Julia validate_solution / SolverTolerances",
    )
end

function _final_quality_report(
    provisional::NEGFSolution,
    outer_converged::Bool;
    candidate = nothing,
    audit_failure = nothing,
)
    strict_report = validate_solution(provisional; candidate, audit_failure)
    strict = outer_converged && provisional.scba.converged && strict_report.passed
    strict && return (; report = strict_report, converged = true, approximate = false)
    if outer_converged &&
       scba_accepted(provisional.scba) &&
       provisional.options.convergence.mode === :adaptive_working &&
       provisional.options.convergence.diagnostic_quality.enabled
        candidate = NEGFSolution(
            provisional.problem,
            _approximate_options(provisional.options),
            provisional.Uᴴ,
            provisional.n,
            provisional.scba,
            provisional.outer_history,
            provisional.observables,
            provisional.report,
            false,
            provisional.status,
        )
        approximate_report = _validate_solution_metrics(
            candidate,
            strict_report.metrics;
            accept_approximate_scba = true,
        )
        if approximate_report.passed
            return (; report = approximate_report, converged = false, approximate = true)
        end
    end
    return (; report = strict_report, converged = false, approximate = false)
end

"""Propagate final scientific acceptance to the sweep, report and queue warnings."""
function _final_validation_warning(provisional::NEGFSolution, assessment)
    assessment.converged && return nothing
    # Failed/incomplete states remain failed. This warning only makes the
    # physical reason available without a binary checkpoint or solver rerun.
    approximate = assessment.approximate
    options = approximate ? _approximate_options(provisional.options) : provisional.options
    limits = _solution_validation_limits(options.tolerances, assessment.report.metrics)
    return Dict{String,Any}(
        "code" => approximate ? "SOLUTION_APPROXIMATE_ACCEPTED" : "FINAL_VALIDATION_FAILED",
        "scope" => "point",
        "message" =>
            approximate ?
            "Final state passed the declared approximate band; strict convergence is not certified." :
            "Final scientific acceptance failed: " * (
                isempty(assessment.report.messages) ?
                "Outer Poisson fixed point did not reach its declared stopping condition." :
                join(assessment.report.messages, "; ")
            ),
        "thresholds" => Dict(String(name)=>value for (name, value) in limits),
        "metrics" => Dict(
            String(name)=>isfinite(value) ? value : nothing for
            (name, value) in assessment.report.metrics
        ),
    )
end

"""Independent hot-LO power check: phonon event counts versus stored electron collision energy."""
function _lo_balance_diagnostics(problem::NEGFProblem, scba::SCBAResult)
    kernel = problem.kernels.K[:LO]
    scale = problem.kernels.qᴷ[:LO]
    emission = _static_contraction(
        kernel,
        apply_energy_shift(problem.W₋ᴸᴼ, scba.green.Gᵍ),
        problem.grids.wᵏ,
        scale,
    )
    absorption = _static_contraction(
        kernel,
        apply_energy_shift(problem.W₊ᴸᴼ, scba.green.Gᵍ),
        problem.grids.wᵏ,
        scale,
    )
    state = lo_kinetic_state(problem, scba.green, emission, absorption)
    electron_gain = Float64(
        ustrip(u"W/m^2", _collision_power(problem, scba.scattering[:LO], scba.green)),
    )
    residual =
        abs(electron_gain+state.electron_power_W_per_m2) /
        max(abs(electron_gain), abs(state.electron_power_W_per_m2), 1e-30)
    return (;
        state,
        electron_gain_W_per_m2 = electron_gain,
        collision_balance_residual = residual,
        rates_source = :fresh_unmixed_literal_map,
        collision_source = :published_stored_family,
    )
end

function _publish_model_diagnostics!(observables, problem::NEGFProblem, scba::SCBAResult, U)
    observables[:physical_models] = _restart_contract_value(problem.models)
    observables[:effective_seed] = Dict{String,Any}(
        String(name) => (value isa Symbol ? String(value) : value) for
        (name, value) in pairs(effective_seed_parameters(problem))
    )
    measured = measured_representation_coverage(problem, scba, U)
    estimate = measured.estimate
    observables[:representation_coverage] = Dict{String,Any}(
        "source"=>"measured_state",
        "discretization_id"=>estimate.discretization_id,
        "basis_id"=>estimate.basis_id,
        "hartree_id"=>estimate.hartree_id,
        "hamiltonian_covered"=>estimate.covered,
        "energy_unit"=>"E0",
        "energy_reference_eV"=>_electronvolts(problem.physical.E_ref),
        "E0_eV"=>problem.scales.E₀_eV,
        "energy_min"=>estimate.energy_min,
        "energy_max"=>estimate.energy_max,
        "hamiltonian_min"=>estimate.h_min,
        "hamiltonian_max"=>estimate.h_max,
        "required_min"=>estimate.required_min,
        "required_max"=>estimate.required_max,
        "spectral_matrix_deficit"=>measured.spectral_matrix_deficit,
        "energy_step"=>measured.energy_step,
        "minimum_resolved_broadening"=>measured.minimum_resolved_broadening,
        "tails"=>Dict(
            String(name)=>getfield(measured.tails, name) for name in keys(measured.tails)
        ),
    )
    if :electron_electron in problem.kernels.enabled
        diagnostics = electron_electron_collision_diagnostics(
            problem,
            scba.green,
            scba.scattering[:electron_electron],
        )
        observables[:electron_electron] =
            Dict(String(name)=>value for (name, value) in pairs(diagnostics))
    end
    if problem.scattering.LO
        diagnostics = try
            _lo_balance_diagnostics(problem, scba)
        catch error
            error isa DomainError || rethrow()
            observables[:lo_population] =
                Dict("available"=>false, "reason"=>sprint(showerror, error))
            nothing
        end
        if diagnostics !== nothing
            state = diagnostics.state
            data = Dict{String,Any}(
                String(name)=>getfield(state, name) for name in fieldnames(LOKineticState)
            )
            data["electron_gain_W_per_m2"] = diagnostics.electron_gain_W_per_m2
            data["collision_balance_residual"] = diagnostics.collision_balance_residual
            data["rates_source"] = String(diagnostics.rates_source)
            data["collision_source"] = String(diagnostics.collision_source)
            observables[:lo_population] = data
        end
    end
    return observables
end
