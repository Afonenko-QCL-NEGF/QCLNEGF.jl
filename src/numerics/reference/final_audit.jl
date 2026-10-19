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
fresh_fixed_point_metrics(problem::NEGFProblem, scba::SCBAResult,
                          U::AbstractVector{<:Real}, candidate::FixedPointAuditCandidate) =
    _fresh_fixed_point_measurements(problem, scba, U, candidate).metrics

function _fresh_fixed_point_measurements(
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
    sigma = Dict{Symbol,Any}(
        :embedding_total=>_selfenergy_residual_evidence(scba.embedding, candidate.embedding),
        :embedding_plus=>_selfenergy_residual_evidence(scba.embedding_plus, candidate.plus),
        :embedding_minus=>_selfenergy_residual_evidence(scba.embedding_minus, candidate.minus),
    )
    r_sigma = maximum((sigma[:embedding_total].ratio, sigma[:embedding_plus].ratio,
                       sigma[:embedding_minus].ratio))
    for name in problem.kernels.enabled
        sigma[name] = _selfenergy_residual_evidence(scba.scattering[name], candidate.scattering[name])
        r_sigma = max(r_sigma, sigma[name].ratio)
    end
    dyson = _dyson_residual_evidence(problem.grids.ε, h, candidate_total.Σᴿ, scba.green.Gᴿ)
    spectral = _spectral_residual_evidence(scba.green, candidate_total)
    keldysh = _keldysh_residual_evidence(scba.green, candidate_total)
    metrics = Dict{Symbol,Float64}(
        :r_D=>_dyson_residual(problem.grids.ε, h, total.Σᴿ, scba.green.Gᴿ),
        :r_D_candidate=>dyson.ratio,
        :r_A_candidate=>spectral.ratio,
        :r_A=>_spectral_residual(scba.green, total),
        :r_K=>keldysh.ratio,
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
    return (; metrics, evidence=(; dyson, spectral, keldysh, sigma))
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
        "final_fixed_point_audit"=>get(solution.observables, :final_fixed_point_audit, nothing),
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

# SCI05 internal template stage: no execution identity is invented here. Runner
# binds its six actual state_identity fields before validating/serializing the
# portable v1 record. This template alone is not a schema-valid final artifact.
const _FINAL_AUDIT_SCHEMA = "qcl-negf-final-fixed-point-audit-v1"
const _FINAL_AUDIT_CHANNELS = (:LO, :acoustic, :impurity, :IFR, :alloy, :electron_electron)
const _FINAL_AUDIT_FLOOR = Float64(1e-14)

function _audit_measure(value::Real; invalid_reason=nothing)
    x = Float64(value)
    tag = isnan(x) ? "nan" : x == Inf ? "+inf" : x == -Inf ? "-inf" : nothing
    return Dict{String,Any}(
        "status"=>tag !== nothing || invalid_reason !== nothing ? "error" : "measured",
        "value"=>tag === nothing ? (iszero(x) ? 0.0 : x) : nothing,
        "nonfinite"=>tag,
        "reason"=>tag === nothing ? invalid_reason : "nonfinite_measurement",
    )
end

_audit_absent(reason="map_unavailable"; status="unavailable") = Dict{String,Any}(
    "status"=>status, "value"=>nothing, "nonfinite"=>nothing, "reason"=>reason,
)
_audit_undefined() = _audit_absent("invalid_prerequisite"; status="error")

function _audit_detail(text::AbstractString)
    output = IOBuffer()
    count = 0
    truncated = false
    for c in text
        width = ncodeunits(string(c))
        if count + width > 512
            truncated = true
            break
        end
        print(output, c)
        count += width
    end
    return String(take!(output)), truncated
end

function _audit_norm_record(raw)
    record = Dict{String,Any}("floor"=>_FINAL_AUDIT_FLOOR)
    for name in (:numerator_F, :scale_F, :denominator, :ratio)
        record[String(name)] = raw === nothing ? _audit_absent() : _audit_measure(getfield(raw, name))
    end
    return record
end

# D2 status clarification: preserve the old finite quotient even when its
# normalization overflowed. Only actual nonfinite fields receive a tag.
_audit_dyson_ratio_measure(raw, value=raw.ratio) = _audit_measure(value;
    invalid_reason=!isfinite(raw.denominator) ? "nonfinite_measurement" : nothing)

function _audit_witness_record(raw)
    names = (:numerator_F, :D_F, :GR_F, :identity_F, :denominator, :ratio)
    record = Dict{String,Any}(
        "status"=>raw === nothing ? "unavailable" : "measured",
        "e_index"=>raw === nothing ? nothing : raw.e_index,
        "k_index"=>raw === nothing ? nothing : raw.k_index,
    )
    for name in names
        measurement = raw === nothing ? _audit_absent() : name === :ratio ?
                      _audit_dyson_ratio_measure(raw) : _audit_measure(getfield(raw, name))
        record[String(name)] = measurement
        measurement["status"] == "error" && (record["status"] = "error")
    end
    return record
end

function _final_audit_context(solution)
    problem = solution.problem
    producer = get(solution.observables, :acceptance_measurement_context, Dict())
    algorithms = get(get(solution.observables, :restart_contract, Dict()), "algorithms", Dict())
    # Every nullable selector is supplied by its actual owner, never AlgorithmOptions defaults.
    context = Dict{String,Any}(
        "backend"=>get(producer, "backend", nothing),
        "energy_shift"=>get(algorithms, "energy_shift", nothing),
        "hilbert"=>get(algorithms, "hilbert", nothing),
        "contraction"=>get(algorithms, "contraction", nothing),
        "kernel_build"=>get(algorithms, "kernel_build", nothing),
        "embedding"=>get(algorithms, "embedding", nothing),
        "electron_electron_mode"=>String(problem.models.electron_electron.mode),
        "embedding_periods"=>get(algorithms, "embedding_periods", nothing),
        "enabled_scattering"=>String.(_final_audit_channels(problem)),
        "E0_eV"=>problem.scales.E₀_eV,
        "L0_m"=>problem.scales.L₀_m,
        "lambda_P"=>problem.scales.λ_P,
        "map_kernels_shared_with_solver"=>get(producer, "map_kernels_shared_with_solver", nothing),
        "independent_physics_oracle"=>false,
        "state_alignment"=>get(producer, "final_state_alignment", "unknown"),
        "alignment_reason"=>get(producer, "final_state_alignment", "unknown") == "unknown" ?
                            "producer_precondition_unverified" : nothing,
        "source_owner_role"=>"physics.full",
        "source_paths"=>Dict("U"=>"/state_dimensionless/UH", "n"=>"/state_dimensionless/density",
                             "GR"=>"/state_dimensionless/GR", "GL"=>"/state_dimensionless/GL",
                             "GG"=>"/state_dimensionless/GG", "A"=>"/state_dimensionless/A"),
    )
    for name in ("E0_eV", "L0_m", "lambda_P")
        value = context[name]
        (!isfinite(value) || value <= 0) && (context[name] = nothing)
    end
    nullable = ("backend", "energy_shift", "hilbert", "contraction", "kernel_build", "embedding",
                "electron_electron_mode", "embedding_periods", "E0_eV", "L0_m", "lambda_P",
                "map_kernels_shared_with_solver")
    unknown = String[name for name in nullable if context[name] === nothing]
    context["unknown_fields"] = unknown
    context["status"] = isempty(unknown) ? "known" : "unknown"
    context["reason"] = isempty(unknown) ? nothing : "missing_executed_context"
    return context
end

_final_audit_channels(problem) = [name for name in _FINAL_AUDIT_CHANNELS if name in problem.kernels.enabled]

function _audit_history(solution)
    inner, outer = solution.scba.history, solution.outer_history
    previous = length(inner) > 1 ? _audit_measure(inner[end-1].λ) :
               isempty(inner) ? _audit_absent("history_unavailable") :
               _audit_absent("no_previous_inner_iteration"; status="not_applicable")
    return Dict{String,Any}(
        "inner_last"=>isempty(inner) ? nothing : inner[end].ν,
        "inner_previous"=>length(inner) > 1 ? inner[end-1].ν : nothing,
        "outer_last"=>isempty(outer) ? nothing : outer[end].μ,
        "previous_lambda"=>previous,
    )
end

function _audit_charge_record(solution, metrics, available, history)
    problem = solution.problem
    if !available
        # The unavailable path did not measure stored Keldysh numbers. Do not
        # allocate another total/Graw field or repeat the map to fill that gap.
        target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
        record = Dict{String,Any}(name=>_audit_absent() for name in (
            "stored_raw_sheet", "normalized_sheet", "candidate_raw_sheet", "stored_lambda",
            "candidate_lambda", "stored_constraint", "candidate_constraint", "normalized_residual",
        ))
        record["target_sheet"] = _audit_measure(target;
            invalid_reason=target <= 0 ? "nonpositive_target" : nothing)
        previous = history["previous_lambda"]
        record["lambda_change"] = previous["status"] in ("unavailable", "not_applicable") ?
                                  copy(previous) : _audit_absent()
        return record
    end
    target, raw, normalized = metrics[:target_charge], metrics[:raw_charge], metrics[:normalized_charge]
    valid_target = isfinite(target) && target > 0
    valid_raw = isfinite(raw) && raw > 0
    record = Dict{String,Any}(
        "target_sheet"=>_audit_measure(target; invalid_reason=target <= 0 ? "nonpositive_target" : nothing),
        "stored_raw_sheet"=>_audit_measure(raw; invalid_reason=raw <= 0 ? "nonpositive_stored_raw" : nothing),
        "normalized_sheet"=>_audit_measure(normalized),
        "stored_lambda"=>valid_target && valid_raw ? _audit_measure(raw/target) : _audit_undefined(),
        "stored_constraint"=>valid_target && valid_raw ? _audit_measure(abs(raw/target-1)) : _audit_undefined(),
        "normalized_residual"=>valid_target && isfinite(normalized) ?
                               _audit_measure(abs(normalized-target)/target) : _audit_undefined(),
    )
    if available
        candidate_raw = metrics[:candidate_raw_charge]
        record["candidate_raw_sheet"] = _audit_measure(candidate_raw)
        valid_candidate = valid_target && isfinite(candidate_raw)
        record["candidate_lambda"] = valid_candidate ? _audit_measure(metrics[:candidate_lambda]) : _audit_undefined()
        record["candidate_constraint"] = valid_candidate ? _audit_measure(metrics[:r_lambda_candidate]) : _audit_undefined()
    else
        for name in ("candidate_raw_sheet", "candidate_lambda", "candidate_constraint")
            record[name] = _audit_absent()
        end
    end
    previous = history["previous_lambda"]
    record["lambda_change"] = previous["status"] == "measured" ?
        (valid_target && valid_raw ? _audit_measure(raw/target-previous["value"]) : _audit_undefined()) :
        previous["status"] == "error" ? _audit_undefined() : copy(previous)
    return record
end

function _audit_has_error(value)
    value isa AbstractDict || value isa AbstractVector || return false
    value isa AbstractDict && get(value, "status", nothing) == "error" && return true
    entries = value isa AbstractDict ? values(value) : value
    return any(_audit_has_error, entries)
end

function _final_fixed_point_record(solution, evidence, metrics, poisson; audit_failure=nothing)
    available = evidence !== nothing
    NE, Nk, Nb, _ = size(solution.scba.green.Gᴿ)
    detail, detail_truncated = _audit_detail(audit_failure === nothing ?
        (available ? "" : String(solution.scba.status)) : String(audit_failure))
    dyson = Dict{String,Any}()
    for name in (:ratio, :global_abs_F, :max_abs_F)
        dyson[String(name)] = !available ? _audit_absent() : name === :ratio ?
                             _audit_dyson_ratio_measure(evidence.dyson.ratio_witness, evidence.dyson.ratio) :
                             _audit_measure(getfield(evidence.dyson, name))
    end
    for name in (:ratio_witness, :absolute_witness)
        dyson[String(name)] = _audit_witness_record(available ? getfield(evidence.dyson, name) : nothing)
    end
    rows = Dict{String,Any}[]
    winner, best = nothing, 0.0
    for family in (_final_audit_channels(solution.problem)..., :embedding_total, :embedding_plus, :embedding_minus)
        for (index, component) in enumerate(("R", "L", "G"))
            row = Dict{String,Any}("family"=>String(family), "component"=>component)
            raw = available ? evidence.sigma[family].rows[index] : nothing
            for name in (:numerator_F, :scale_F, :denominator, :ratio)
                row[String(name)] = available ? _audit_measure(getfield(raw, name)) : _audit_absent()
            end
            if available && (winner === nothing || raw.ratio > best)
                winner = Dict("family"=>String(family), "component"=>component)
            end
            available && (best = max(best, raw.ratio))
            push!(rows, row)
        end
    end
    any(_audit_has_error, rows) && (winner = nothing)
    selfenergy = Dict{String,Any}(
        "rows"=>rows, "winner"=>winner, "floor"=>_FINAL_AUDIT_FLOOR,
        "max_ratio"=>available ? _audit_measure(metrics[:r_Σ]) : _audit_absent(),
    )
    history = _audit_history(solution)
    charge = _audit_charge_record(solution, metrics, available, history)
    poisson_record = Dict{String,Any}(String(name)=>_audit_measure(value) for (name, value) in pairs(poisson))
    poisson_record["rU_max"] = _audit_measure(metrics[:r_U])
    if isfinite(poisson.target_sheet) && poisson.target_sheet <= 0
        poisson_record["target_sheet"] = _audit_measure(poisson.target_sheet; invalid_reason="nonpositive_target")
    end
    if !isfinite(poisson.target_sheet) || poisson.target_sheet <= 0
        poisson_record["neutrality_ratio"] = _audit_undefined()
    end
    spectral = _audit_norm_record(available ? evidence.spectral : nothing)
    keldysh = _audit_norm_record(available ? evidence.keldysh : nothing)
    # Optional history is excluded from the mandatory measurement status.
    mandatory_charge = Dict(name=>value for (name,value) in charge if name != "lambda_change")
    failed = any(_audit_has_error, (dyson, spectral, keldysh, selfenergy, mandatory_charge, poisson_record))
    return Dict{String,Any}(
        "schema"=>_FINAL_AUDIT_SCHEMA,
        "record_status"=>!available ? "unavailable" : failed ? "partial_error" : "measured",
        "reason"=>!available ? "map_unavailable" : failed ? "measurement_error" : nothing,
        "detail"=>detail, "detail_truncated"=>detail_truncated,
        "shape"=>Dict("NE"=>NE, "Nk"=>Nk, "Nb"=>Nb),
        "units"=>Dict("native"=>"1", "energy_scale"=>"eV", "length_scale"=>"m", "index"=>"1-based"),
        "context"=>_final_audit_context(solution), "history"=>history,
        "dyson"=>dyson, "spectral"=>spectral, "keldysh"=>keldysh,
        "selfenergy"=>selfenergy, "charge"=>charge, "poisson"=>poisson_record,
    )
end
