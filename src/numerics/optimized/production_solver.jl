"""
    solve_scba_production(problem, Uᴴ; ...)

Fixed-Hartree SCBA with exact sparse shifts, chunked BLAS contractions,
zero-padded FFT Hilbert transforms, and threaded independent Dyson/Keldysh
blocks.  It returns the same `SCBAResult` contract as the educational direct
backend.  All mechanism candidates are still built from one common Green
state before mixing (Jacobi semantics).

`iteration_callback`, when supplied, is invoked synchronously with a borrowed
view of the accepted arrays before in-place mixing.  The callback may inspect
or serialize that state during the call, but must make its own copy if it needs
to retain the object after returning. `consume_initial=true` transfers ownership
of a private restart's arrays and permits in-place mixing; callers must not reuse
that initial state. The default preserves caller-owned input arrays.

See [the SCBA fixed point](@ref theory-scba),
[Microscopic kernels](@ref theory-kernels), and
[Production backend](@ref theory-production).
"""
function solve_scba_production(
    problem::NEGFProblem,
    Uᴴ::AbstractVector{<:Real};
    options::SolverOptions = baseline_options(),
    production_options::ProductionOptions = ProductionOptions(),
    cache::Union{Nothing,ProductionCache} = nothing,
    initial::Union{Nothing,SCBAResult} = nothing,
    consume_initial::Bool = false,
    iteration_callback::Union{Nothing,Function} = nothing,
    history_observer::Union{Nothing,Function} = nothing,
    resume_iterations::Bool = false,
    restart_options::SolverOptions = options,
)
    _check_options(options)
    _check_production_options(production_options)
    owns_scba_stage = _begin_solver_stage(
        production_options,
        :scba;
        label = "fixed Hartree field",
        iteration = resume_iterations && initial !== nothing && !isempty(initial.history) ? last(initial.history).ν : 0,
        total = options.max_scba,
    )
    cache === nothing &&
        (cache = build_production_cache(problem; options = production_options))
    cache = _validated_production_cache(cache, problem, production_options)
    h = project_hamiltonians(problem, Uᴴ)
    cavity_plan =
        production_options.algorithms.embedding === :finite_chain ?
        _cavity_production_plan(problem, h, production_options) : nothing
    seed_mu_scaled = NaN
    seed_number_ratio = NaN
    n = problem.numerical
    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    field_bytes = sizeof(ComplexF64)*prod(shape)
    history =
        resume_iterations && initial !== nothing ? copy(initial.history) : SCBAIteration[]
    completed_iterations = isempty(history) ? 0 : last(history).ν
    if !isempty(history)
        for previous in history
            marker = previous.physical_markers
            if marker !== nothing && isfinite(marker.seed_mu_eV)
                seed_mu_scaled = marker.seed_mu_eV/problem.scales.E₀_eV
                seed_number_ratio = marker.seed_number_ratio
                break
            end
        end
    end
    completed_iterations <= options.max_scba ||
        throw(ArgumentError("checkpoint exceeds the configured SCBA budget"))
    if resume_iterations && initial !== nothing
        _check_exact_restart_contract(
            initial.restart_contract,
            restart_options,
            production_options.algorithms,
        )
    end
    residual_jobs = cld(n.N_E * n.N_k, production_options.residual_chunk)
    residual_workers =
        production_options.parallel_backend === :blas ? 1 :
        _production_thread_worker_count(production_options, residual_jobs)
    residual_workspace = ProductionResidualWorkspace(
        problem;
        chunk_size = production_options.residual_chunk,
        worker_count = residual_workers,
    )
    if production_options.parallel_backend === :threads &&
       _production_thread_worker_count(
           production_options,
           cld(n.N_E, production_options.energy_chunk),
       ) > 1 &&
       BLAS.get_num_threads() != 1
        @warn "threaded production contractions should use one BLAS thread to avoid nested oversubscription" julia_threads=Base.Threads.nthreads(
            :default,
        ) production_workers=_production_thread_worker_count(
            production_options,
            cld(n.N_E, production_options.energy_chunk),
        ) blas_threads=BLAS.get_num_threads()
    end
    had_initial = initial !== nothing
    initialization_started = time_ns()
    initialization_phase = _production_phase_begin(
        production_options,
        :initialization,
        completed_iterations;
        task_width = _production_thread_worker_count(production_options, n.N_E*n.N_k),
        memory_burst_bytes = initial === nothing ?
                             (24+3length(problem.kernels.enabled))*field_bytes :
                             consume_initial ? 0 :
                             (9+3length(problem.kernels.enabled))*field_bytes,
    )
    initialization_status = :failed
    local green, scattering, embedding, plus, minus
    try
        if initial === nothing
            green, seed_mu_scaled =
                _seed_green_production(problem, h, initialization_phase.options)
            seed_number_ratio =
                _number_functional_production(
                    green.Gˡ,
                    problem.grids,
                    problem.physical.g_s,
                    initialization_phase.options,
                ) / _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
            scattering = _scattering_candidate_production(
                problem,
                green,
                cache,
                initialization_phase.options,
            )
            embedding, plus, minus = _embedding_self_energy_selected(
                problem,
                green,
                cache,
                initialization_phase.options;
                scattering,
                h,
                cavity_plan,
            )
        else
            Set(keys(initial.scattering)) == Set(problem.kernels.enabled) ||
                throw(ArgumentError("initial SCBA mechanisms do not match problem"))
            all(
                size(getfield(family, component)) == shape for
                family in values(initial.scattering) for component in (:Σᴿ, :Σˡ, :Σᵍ)
            ) || throw(DimensionMismatch("initial scattering arrays have wrong shape"))
            all(
                size(getfield(family, component)) == shape for family in
                (initial.embedding, initial.embedding_plus, initial.embedding_minus) for
                component in (:Σᴿ, :Σˡ, :Σᵍ)
            ) || throw(DimensionMismatch("initial embedding arrays have wrong shape"))
            all(
                size(getfield(initial.green, component)) == shape for
                component in (:Gᴿ, :Gˡ, :Gᵍ, :A)
            ) || throw(DimensionMismatch("initial Green arrays have wrong shape"))
            size(initial.green.condition_number) == (n.N_E, n.N_k) &&
            size(initial.green.dyson_scale) == (n.N_E, n.N_k) ||
                throw(DimensionMismatch("initial Green diagnostics have wrong shape"))
            green = initial.green
            # Freshly loaded private restart buffers may transfer ownership.
            # Public/warm inputs retain the historical non-mutating contract.
            scattering =
                consume_initial ? initial.scattering : _copy_mechanisms(initial.scattering)
            embedding =
                consume_initial ? initial.embedding : _copy_family(initial.embedding)
            plus =
                consume_initial ? initial.embedding_plus :
                _copy_family(initial.embedding_plus)
            minus =
                consume_initial ? initial.embedding_minus :
                _copy_family(initial.embedding_minus)
        end
        initialization_status = :completed
    finally
        _production_phase_end(
            production_options,
            initialization_phase;
            status = initialization_status,
        )
    end
    if production_options.phase_timing && owns_scba_stage
        _update_solver_stage(
            production_options,
            :scba;
            iteration = completed_iterations,
            total = options.max_scba,
            metrics = SolverMetric[
                SolverMetric(
                    :t_initialization,
                    (time_ns() - initialization_started) * 1e-9,
                    "s",
                ),
                SolverMetric(:seeded, initial === nothing),
            ],
            message = "initial state ready",
        )
    end
    convergence_streak = _trailing_pass_count(
        history,
        row ->
            scba_convergence_assessment(row, options.tolerances, options.convergence).passed,
    )
    diagnostic_quality_streak =
        options.convergence.mode === :adaptive_working &&
        options.convergence.diagnostic_quality.enabled ?
        _trailing_pass_count(
            history,
            row -> scba_diagnostic_quality_assessment(
                row,
                options.tolerances,
                options.convergence,
            ).passed,
        ) : 0
    if completed_iterations > 0
        # Checkpoints are aligned with the current state and emitted before
        # terminal decisions. Replay that decision before another mixing step.
        restored_strict =
            convergence_streak >= options.convergence.required_consecutive_scba_passes
        restored_approximate =
            !restored_strict && scba_approximate_acceptance(history, options)
        restored_stop =
            _scba_stopping_reason(history, options.tolerances, options.convergence)
        restored_stagnated =
            !restored_strict && !restored_approximate && restored_stop !== nothing
        if restored_strict || restored_approximate || restored_stagnated
            restored_total = _total_family_production(
                scattering,
                embedding,
                problem.kernels.enabled,
                production_options,
            )
            _ensure_scba_physical_markers!(
                history,
                problem,
                green,
                restored_total,
                production_options;
                seed_mu_scaled,
                seed_number_ratio,
                terminal = true,
            )
            restored_status =
                restored_strict ? :converged :
                restored_approximate ? :approximate : restored_stop
            restored_quality =
                restored_strict ? :strictly_converged :
                restored_approximate ? :approximate_fixed_point : :unresolved
            _end_solver_stage(
                production_options,
                :scba,
                owns_scba_stage;
                status = restored_stagnated ? :incomplete : :completed,
                iteration = completed_iterations,
                total = options.max_scba,
                message = "restored SCBA termination decision",
            )
            return SCBAResult(
                green,
                scattering,
                embedding,
                plus,
                minus,
                history,
                restored_strict,
                restored_status,
                restored_quality,
                initial.restart_contract,
                initial.mixer_state,
            )
        end
    end
    if completed_iterations == options.max_scba
        initial === nothing && error("missing SCBA restart state")
        restored_total = _total_family_production(
            scattering,
            embedding,
            problem.kernels.enabled,
            production_options,
        )
        _ensure_scba_physical_markers!(
            history,
            problem,
            green,
            restored_total,
            production_options;
            seed_mu_scaled,
            seed_number_ratio,
            terminal = true,
        )
        _end_solver_stage(
            production_options,
            :scba,
            owns_scba_stage;
            status = :incomplete,
            iteration = completed_iterations,
            total = options.max_scba,
            message = "checkpoint has exhausted the SCBA budget",
        )
        decision = research_continuation_decision(
            history,
            options;
            usable = initial.quality !== :invalid,
        )
        transition = decision.action === :continue_poisson
        return SCBAResult(
            green,
            scattering,
            embedding,
            plus,
            minus,
            history,
            false,
            transition ? :research_continue : :max_iterations,
            transition ? :unresolved :
            initial.quality === :strictly_converged ? :unresolved : initial.quality,
            initial.restart_contract,
            initial.mixer_state,
        )
    end
    witness_retention = _SCBAWitnessRetention(history)
    anderson_workspace = if resume_iterations && initial !== nothing
        initial.mixer_state.method == production_options.algorithms.mixing || throw(
            ArgumentError("checkpoint mixer method differs from the requested algorithm"),
        )
        _AndersonWorkspace(
            consume_initial ? initial.mixer_state.states :
            deepcopy(initial.mixer_state.states),
            consume_initial ? initial.mixer_state.residuals :
            deepcopy(initial.mixer_state.residuals),
        )
    else
        _AndersonWorkspace()
    end
    mixer_state() = SCBAMixerState(
        production_options.algorithms.mixing,
        copy(anderson_workspace.states),
        copy(anderson_workspace.residuals),
    )
    restart_contract =
        _solver_restart_contract(restart_options, production_options.algorithms)
    # Only transferred/copied arrays are needed below. Do not retain the old
    # aggregate Green state throughout all subsequent iterations.
    initial = nothing
    if completed_iterations > 0
        resume_phase = _production_phase_begin(
            production_options,
            :restart_scattering,
            completed_iterations;
            task_width = _production_worker_count(
                production_options,
                cld(n.N_E, production_options.energy_chunk),
            ),
            memory_burst_bytes = (3length(problem.kernels.enabled)+6)*field_bytes,
        )
        next_scattering =
            _scattering_candidate_production(problem, green, cache, resume_phase.options)
        _production_phase_end(production_options, resume_phase)
        resume_phase = _production_phase_begin(
            production_options,
            :restart_embedding,
            completed_iterations;
            task_width = _production_worker_count(production_options, n.N_E*n.N_k),
            memory_burst_bytes = 12field_bytes,
        )
        next_embedding, next_plus, next_minus = _embedding_self_energy_selected(
            problem,
            green,
            cache,
            resume_phase.options;
            scattering = next_scattering,
            h,
            cavity_plan,
        )
        # The summed embedding is not used in this reconstruction.
        next_embedding = nothing
        _production_phase_end(production_options, resume_phase)
        history_growth =
            production_options.algorithms.mixing === :anderson &&
            length(anderson_workspace.states)+length(anderson_workspace.free_states) <
            production_options.algorithms.anderson_history_depth
        resume_phase = _production_phase_begin(
            production_options,
            :restart_mixing,
            completed_iterations;
            task_width = _production_thread_worker_count(production_options, prod(shape)),
            memory_burst_bytes = (
                history_growth ? 6(length(problem.kernels.enabled)+2)*field_bytes : 0
            ) +
                                 128n.N_b^2*sizeof(ComplexF64)*_production_thread_worker_count(
                production_options,
                prod(shape),
            ),
        )
        if production_options.algorithms.mixing === :anderson
            _anderson_mix!(
                anderson_workspace,
                _ordered_iteration_families(
                    scattering,
                    plus,
                    minus,
                    problem.kernels.enabled,
                ),
                _ordered_iteration_families(
                    next_scattering,
                    next_plus,
                    next_minus,
                    problem.kernels.enabled,
                ),
                production_options.algorithms,
                options.α_Σ,
                resume_phase.options,
            )
        else
            _mix_mechanisms_production!(
                scattering,
                next_scattering,
                problem.kernels.enabled,
                options.α_Σ,
                resume_phase.options,
            )
            _mix_family_production!(plus, next_plus, options.α_Σ, resume_phase.options)
            _mix_family_production!(minus, next_minus, options.α_Σ, resume_phase.options)
        end
        _embedding_sum_production!(embedding, plus, minus, resume_phase.options)
        _production_phase_end(production_options, resume_phase)
        next_scattering = next_plus = next_minus = nothing
    end
    Jprevious = isempty(history) ? nothing : last(history).J
    populations_previous =
        !had_initial ? nothing :
        real.(diag(_sheet_density_matrix_bar(problem, green.Gˡ))) ./
        _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    for ν = (completed_iterations+1):options.max_scba
        iteration_started = time_ns()
        dyson_started = iteration_started
        phase = _production_phase_begin(
            production_options,
            :dyson,
            ν;
            task_width = _production_worker_count(production_options, n.N_E*n.N_k),
            work_units = n.N_E*n.N_k,
            memory_burst_bytes = 8field_bytes,
        )
        local green_current, λ, total
        try
            green_current, λ, total =
                _current_green_production(problem, h, scattering, embedding, phase.options)
            green = green_current
        catch error
            _production_phase_end(production_options, phase; status = :failed)
            _end_solver_stage(
                production_options,
                :scba,
                owns_scba_stage;
                status = :failed,
                iteration = ν,
                total = options.max_scba,
                message = "Dyson/Keldysh failure",
            )
            @warn "SCBA stopped during Dyson/Keldysh at iteration $ν" exception=(
                error,
                catch_backtrace(),
            )
            # No Green state aligned with the requested Hartree field exists;
            # returning the warm state would make a corrupt restart possible.
            rethrow()
        end
        _production_phase_end(production_options, phase)
        dyson_seconds = (time_ns() - dyson_started) * 1e-9
        phase = _production_phase_begin(
            production_options,
            :candidate_scattering,
            ν;
            task_width = _production_worker_count(
                production_options,
                cld(n.N_E, production_options.energy_chunk),
            ),
            work_units = n.N_E*n.N_k*n.N_b^2,
            memory_burst_bytes = (3length(problem.kernels.enabled)+6)*field_bytes,
        )
        candidate_started = time_ns()
        local candidate_sc, candidate_embedding, candidate_plus, candidate_minus, rround
        local candidate_scattering_seconds, candidate_embedding_seconds
        try
            candidate_sc, rround = _scattering_candidate_production(
                problem,
                green,
                cache,
                phase.options;
                return_roundoff = true,
            )
            _production_phase_end(production_options, phase)
            candidate_scattering_seconds = (time_ns()-candidate_started)*1e-9
            phase = _production_phase_begin(
                production_options,
                :candidate_embedding,
                ν;
                task_width = _production_worker_count(production_options, n.N_E*n.N_k),
                work_units = n.N_E*n.N_k,
                memory_burst_bytes = 12field_bytes,
                workspace_bytes = (
                    8n.N_b^2*sizeof(ComplexF64)+64n.N_b*sizeof(ComplexF64)+(n.N_b+1)*sizeof(
                        LinearAlgebra.BlasInt,
                    )
                )*_production_worker_count(production_options, n.N_E*n.N_k),
            )
            embedding_started = time_ns()
            candidate_embedding, candidate_plus, candidate_minus =
                _embedding_self_energy_selected(
                    problem,
                    green,
                    cache,
                    phase.options;
                    scattering = candidate_sc,
                    h,
                    cavity_plan,
                )
            _production_phase_end(production_options, phase)
            candidate_embedding_seconds = (time_ns()-embedding_started)*1e-9
        catch error
            _production_phase_end(production_options, phase; status = :failed)
            error isa DomainError || rethrow()
            _end_solver_stage(
                production_options,
                :scba,
                owns_scba_stage;
                status = :failed,
                iteration = ν,
                total = options.max_scba,
                message = "candidate construction failure",
            )
            @warn "SCBA candidate failed at iteration $ν" exception=(
                error,
                catch_backtrace(),
            )
            _record_candidate_failure!(history, problem, h, green, total, plus, ν, λ)
            _ensure_scba_physical_markers!(
                history,
                problem,
                green,
                total,
                production_options;
                seed_mu_scaled,
                seed_number_ratio,
                terminal = true,
            )
            history_observer === nothing || history_observer(last(history))
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
            )
        end
        assembly_started = time_ns()
        phase = _production_phase_begin(
            production_options,
            :candidate_assembly,
            ν;
            task_width = _production_thread_worker_count(production_options, prod(shape)),
            work_units = prod(shape),
            memory_burst_bytes = 3field_bytes,
        )
        candidate_total = _total_family_production(
            candidate_sc,
            candidate_embedding,
            problem.kernels.enabled,
            phase.options,
        )
        _production_phase_end(production_options, phase)
        candidate_assembly_seconds = (time_ns()-assembly_started)*1e-9
        candidate_seconds = (time_ns() - candidate_started) * 1e-9

        residual_started = time_ns()
        phase = _production_phase_begin(
            production_options,
            :residuals,
            ν;
            task_width = residual_workers,
            work_units = n.N_E*n.N_k,
            memory_burst_bytes = 0,
        )
        (; rD, rA, rK, rΣ, rλ, rPSD, rcaus) = _production_fixed_point_measurement(
            problem,
            h,
            green,
            total,
            candidate_total,
            scattering,
            embedding,
            plus,
            minus,
            candidate_sc,
            candidate_embedding,
            candidate_plus,
            candidate_minus,
            λ,
            residual_workspace;
            worker_count = phase.task_width,
        )
        _production_phase_end(production_options, phase)
        residual_seconds = (time_ns() - residual_started) * 1e-9

        observables_started = time_ns()
        phase = _production_phase_begin(
            production_options,
            :observables,
            ν;
            work_units = n.N_E*n.N_k,
        )
        (; J, target, populations, rJchange, rpopulation) =
            _production_observable_measurement(
                problem,
                green,
                plus,
                Jprevious,
                populations_previous,
            )
        _production_phase_end(production_options, phase)
        observables_seconds = (time_ns() - observables_started) * 1e-9
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
                _scba_physics_witness(
                    _production_positivity_witness(residual_workspace),
                    green,
                    total,
                    scattering,
                    plus,
                    minus,
                ),
                nothing,
            ),
        )
        _retain_scba_witness_payloads!(history, witness_retention)
        if !all(isfinite, (rD, rA, rK, rΣ, rλ, rPSD, rcaus, rround, J, λ))
            _ensure_scba_physical_markers!(
                history,
                problem,
                green,
                total,
                production_options;
                seed_mu_scaled,
                seed_number_ratio,
                terminal = true,
            )
            history_observer === nothing || history_observer(last(history))
            _end_solver_stage(
                production_options,
                :scba,
                owns_scba_stage;
                status = :failed,
                iteration = ν,
                total = options.max_scba,
                message = "non-finite physics metrics",
            )
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
            )
        end
        assessment = scba_convergence_assessment(
            last(history),
            options.tolerances,
            options.convergence,
        )
        fixed_point = _scba_fixed_point_assessment(
            last(history),
            options.tolerances,
            options.convergence,
        )
        positivity_witness = _production_positivity_witness(residual_workspace)
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
        phase = _production_phase_begin(
            production_options,
            :checkpoint_request,
            ν;
            memory_burst_bytes = 0,
        )
        checkpoint_due =
            iteration_callback !== nothing &&
            _production_checkpoint_requested(production_options, :scba, ν)
        _production_phase_end(production_options, phase)
        phase = _production_phase_begin(production_options, :physics_markers, ν)
        physics_markers_seconds = _ensure_scba_physical_markers!(
            history,
            problem,
            green,
            total,
            production_options;
            seed_mu_scaled,
            seed_number_ratio,
            terminal = decision.terminal || checkpoint_due,
            mixed_scattering = scattering,
            fresh_scattering = candidate_sc,
            mixed_embedding = embedding,
            fresh_embedding = candidate_embedding,
            mixed_plus = plus,
            fresh_plus = candidate_plus,
            mixed_minus = minus,
            fresh_minus = candidate_minus,
            fresh_total = candidate_total,
        )
        _production_phase_end(production_options, phase)
        history_observer === nothing || history_observer(last(history))
        if checkpoint_due
            aligned = SCBAResult(
                green,
                scattering,
                embedding,
                plus,
                minus,
                copy(history),
                false,
                :running,
                quality,
                restart_contract,
                mixer_state(),
            )
            iteration_callback(aligned)
        end
        if decision.terminal
            _emit_scba_progress(
                production_options,
                ν,
                options.max_scba;
                rD,
                rA,
                rK,
                rΣ,
                rλ,
                rPSD,
                fixed_point,
                positivity_witness,
                rcaus,
                rround,
                rJchange,
                rpopulation,
                current = J,
                assessment,
                convergence_streak,
                diagnostic,
                diagnostic_quality_streak,
                quality = decision.quality,
                dyson_seconds,
                candidate_seconds,
                candidate_scattering_seconds,
                candidate_embedding_seconds,
                candidate_assembly_seconds,
                physics_markers_seconds,
                residual_seconds,
                observables_seconds,
                mixing_seconds = 0.0,
                total_seconds = (time_ns()-iteration_started)*1e-9,
            )
            _end_solver_stage(
                production_options,
                :scba,
                owns_scba_stage;
                status = decision.status in (:converged, :approximate, :research_continue) ?
                         :completed : :incomplete,
                iteration = ν,
                total = options.max_scba,
                metrics = SolverMetric[
                    SolverMetric(:J, _solver_current_density(J), "A/cm^2"),
                    SolverMetric(:limiting_gate, String(assessment.limiting_metric)),
                    SolverMetric(:limiting_ratio, assessment.limiting_ratio),
                    SolverMetric(:quality, String(decision.quality)),
                    SolverMetric(:continuation_policy, String(options.convergence.mode)),
                    SolverMetric(:termination_reason, String(decision.status)),
                ],
                message = decision.reason,
            )
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
                mixer_state(),
            )
        end
        mixing_started = time_ns()
        history_growth =
            production_options.algorithms.mixing === :anderson &&
            length(anderson_workspace.states)+length(anderson_workspace.free_states) <
            production_options.algorithms.anderson_history_depth
        mixer_burst =
            (history_growth ? 6(length(problem.kernels.enabled)+2)*field_bytes : 0) +
            128n.N_b^2*sizeof(ComplexF64)*_production_thread_worker_count(
                production_options,
                prod(shape),
            )
        phase = _production_phase_begin(
            production_options,
            :mixing,
            ν;
            task_width = _production_thread_worker_count(production_options, prod(shape)),
            work_units = prod(shape),
            memory_burst_bytes = mixer_burst,
        )
        if production_options.algorithms.mixing === :anderson
            state_families = _ordered_iteration_families(
                scattering,
                plus,
                minus,
                problem.kernels.enabled,
            )
            candidate_families = _ordered_iteration_families(
                candidate_sc,
                candidate_plus,
                candidate_minus,
                problem.kernels.enabled,
            )
            _anderson_mix!(
                anderson_workspace,
                state_families,
                candidate_families,
                production_options.algorithms,
                options.α_Σ,
                phase.options,
            )
        else
            _mix_mechanisms_production!(
                scattering,
                candidate_sc,
                problem.kernels.enabled,
                options.α_Σ,
                phase.options,
            )
            _mix_family_production!(plus, candidate_plus, options.α_Σ, phase.options)
            _mix_family_production!(minus, candidate_minus, options.α_Σ, phase.options)
        end
        _embedding_sum_production!(embedding, plus, minus, phase.options)
        _production_phase_end(production_options, phase)
        mixing_seconds = (time_ns() - mixing_started) * 1e-9
        Jprevious = J
        populations_previous = populations
        total_seconds = (time_ns() - iteration_started) * 1e-9
        _emit_scba_progress(
            production_options,
            ν,
            options.max_scba;
            rD,
            rA,
            rK,
            rΣ,
            rλ,
            rPSD,
            fixed_point,
            positivity_witness,
            rcaus,
            rround,
            rJchange,
            rpopulation,
            current = J,
            assessment,
            convergence_streak,
            diagnostic,
            diagnostic_quality_streak,
            quality,
            dyson_seconds,
            candidate_seconds,
            candidate_scattering_seconds,
            candidate_embedding_seconds,
            candidate_assembly_seconds,
            physics_markers_seconds,
            residual_seconds,
            observables_seconds,
            mixing_seconds,
            total_seconds,
        )
    end
    error("unreachable production SCBA loop exit")
end


"""
    _solve_production(problem; ...)

Fully nested Poisson--SCBA solve using the production inner backend.  Optional
`initial_Uᴴ` and `initial_scba` implement in-memory warm starts across nearby
bias points. An optional neutral state sink receives a state whose
Green functions and self-energies correspond to the stored Hartree field.

See [the closed Poisson–SCBA loop](@ref theory-outer-loop),
[Production backend](@ref theory-production), and the
[state persistence contract](https://github.com/AfonenkoA/QCLNEGFRunner.jl/blob/main/docs/src/user/results.md).
"""
function _solve_production(
    problem::NEGFProblem;
    options::SolverOptions = baseline_options(),
    production_options::ProductionOptions = ProductionOptions(),
    cache::Union{Nothing,ProductionCache} = nothing,
    initial_Uᴴ::Union{Nothing,AbstractVector{<:Real}} = nothing,
    initial_scba::Union{Nothing,SCBAResult} = nothing,
    consume_initial::Bool = false,
    checkpoint_sink::Union{Nothing,Function} = nothing,
    initial_outer_history::Vector{OuterIteration} = OuterIteration[],
    initial_warnings::Vector{Dict{String,Any}} = Dict{String,Any}[],
    resume_scba::Bool = false,
    state_observer::Union{Nothing,Function} = nothing,
    history_observer::Union{Nothing,Function} = nothing,
)
    _check_options(options)
    _check_production_options(production_options)
    cache === nothing &&
        (cache = build_production_cache(problem; options = production_options))
    cache = _validated_production_cache(cache, problem, production_options)
    Nz = problem.numerical.N_z
    Uᴴ = initial_Uᴴ === nothing ? zeros(Float64, Nz) : Float64.(initial_Uᴴ)
    published_U = copy(Uᴴ)
    length(Uᴴ) == Nz || throw(DimensionMismatch("initial Hartree field differs"))
    nprevious = nothing
    previous_scba = initial_scba
    initial_scba = nothing
    Jprevious = nothing
    populations_previous = nothing
    outer_history = copy(initial_outer_history)
    completed_outer = isempty(outer_history) ? 0 : last(outer_history).μ
    completed_outer <= options.max_poisson ||
        throw(ArgumentError("checkpoint exceeds the configured Poisson budget"))
    if !isempty(outer_history)
        nprevious = copy(last(outer_history).density)
        Jprevious = last(outer_history).J
        populations_previous = copy(last(outer_history).populations)
    end
    warnings = deepcopy(initial_warnings)
    inner_policy_history = Dict{String,Any}[]
    restart_contract = _solver_restart_contract(options, production_options.algorithms)
    outer_tolerances = options.tolerances
    convergence_streak = _trailing_pass_count(
        outer_history,
        row ->
            poisson_convergence_assessment(row, outer_tolerances, options.convergence).passed,
    )
    outer_converged =
        convergence_streak >= options.convergence.required_consecutive_poisson_passes &&
        previous_scba !== nothing &&
        previous_scba.converged
    status = outer_converged ? :converged : :max_poisson_iterations
    owns_poisson_stage = _begin_solver_stage(
        production_options,
        :poisson;
        label = "self-consistent Hartree potential",
        iteration = completed_outer,
        total = options.max_poisson,
    )

    function inner_checkpoint_callback(Uaccepted)
        checkpoint_sink === nothing && return nothing
        return function (scba_state::SCBAResult)
            density_phase = _production_phase_begin(
                production_options,
                :checkpoint_density,
                last(scba_state.history).ν;
                stage = :checkpoint,
                work_units = problem.numerical.N_E*problem.numerical.N_k,
                memory_burst_bytes = sizeof(ComplexF64)*problem.numerical.N_b^2 +
                                     sizeof(Float64)*problem.numerical.N_z,
            )
            n̄ = _electron_density_bar(problem, scba_state.green.Gˡ)
            _production_phase_end(production_options, density_phase)
            dummy = ConvergenceReport(false, Dict{Symbol,Float64}(), String[])
            running = NEGFSolution(
                problem,
                options,
                copy(Uaccepted),
                n̄,
                scba_state,
                copy(outer_history),
                Dict{Symbol,Any}(
                    :warnings=>deepcopy(warnings),
                    :restart_contract=>deepcopy(restart_contract),
                ),
                dummy,
                false,
                :running_scba,
            )
            publish_phase = _production_phase_begin(
                production_options,
                :checkpoint_publish,
                last(scba_state.history).ν;
                stage = :checkpoint,
            )
            try
                checkpoint_sink(running)
            finally
                _production_phase_end(production_options, publish_phase)
            end
            return nothing
        end
    end

    for μ = (completed_outer+1):(outer_converged ? completed_outer : options.max_poisson)
        Uaccepted = copy(Uᴴ)
        inner_policy = inner_working_policy(options, outer_history)
        push!(
            inner_policy_history,
            Dict{String,Any}(
                "scope"=>"inner_scba",
                "message"=>"Inner tolerance selected from completed raw outer residuals; final acceptance remains strict.",
                "outer_iteration"=>μ,
                "stage"=>String(inner_policy.stage),
                "forcing_threshold"=>inner_policy.forcing_threshold,
                "outer_raw_residual"=>inner_policy.outer_residual,
                "approximate_enabled"=>inner_policy.options.convergence.diagnostic_quality.enabled,
                "strict_final_required"=>true,
            ),
        )
        scba = solve_scba_production(
            problem,
            Uaccepted;
            options = inner_policy.options,
            # The immutable campaign policy is the exact-restart contract.
            # Its effective inner band is reproducibly derived from the saved
            # completed outer history; it is not a different user input.
            restart_options = options,
            production_options = with_production_options(
                production_options;
                outer_iteration = μ,
            ),
            cache,
            initial = let state = previous_scba
                previous_scba = nothing
                state
            end,
            consume_initial = consume_initial && μ == completed_outer+1,
            iteration_callback = inner_checkpoint_callback(Uaccepted),
            history_observer = history_observer===nothing ? nothing :
                               (row->history_observer(:scba, μ, row, problem)),
            resume_iterations = resume_scba && μ == completed_outer + 1,
        )
        published_U = copy(Uaccepted)
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
        Unew = (1 - options.α_P) .* Uaccepted .+ options.α_P .* candidate
        rU = maximum(abs.(candidate .- Uaccepted))
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
        assessment = poisson_convergence_assessment(
            last(outer_history),
            outer_tolerances,
            options.convergence,
        )
        convergence_streak =
            assessment.passed && scba.converged ? convergence_streak + 1 : 0
        poisson_metrics = SolverMetric[
            SolverMetric(:J, _solver_current_density(J), "A/cm^2"),
            SolverMetric(:r_J, rJ),
            SolverMetric(:r_P, rP),
            SolverMetric(:r_U, rU),
            SolverMetric(:r_n, rn),
            SolverMetric(:limiting_gate, String(assessment.limiting_metric)),
            SolverMetric(:limiting_ratio, assessment.limiting_ratio),
            SolverMetric(:failed_gate_count, length(assessment.failed_metrics)),
            SolverMetric(:convergence_streak, convergence_streak),
        ]
        emitted = _update_solver_stage(
            production_options,
            :poisson;
            iteration = μ,
            total = options.max_poisson,
            metrics = poisson_metrics,
        )
        if !emitted &&
           production_options.progress_every_outer > 0 &&
           (μ == 1 || μ % production_options.progress_every_outer == 0)
            @info "Poisson progress" iteration=μ total=options.max_poisson rP rU rn rJ current_A_per_cm2=_solver_current_density(
                J,
            )
        end

        if checkpoint_sink !== nothing &&
           _production_checkpoint_requested(production_options, :outer, μ)
            dummy = ConvergenceReport(false, Dict{Symbol,Float64}(), String[])
            checkpoint_solution = NEGFSolution(
                problem,
                options,
                Uaccepted,
                n̄,
                scba,
                copy(outer_history),
                Dict{Symbol,Any}(
                    :warnings=>deepcopy(warnings),
                    :restart_contract=>deepcopy(restart_contract),
                ),
                dummy,
                false,
                :running,
            )
            checkpoint_sink(checkpoint_solution)
        end
        if state_observer !== nothing
            state_observer(
                NEGFSolution(
                    problem,
                    options,
                    copy(Uaccepted),
                    n̄,
                    scba,
                    copy(outer_history),
                    Dict{Symbol,Any}(:warnings=>copy(warnings)),
                    ConvergenceReport(false, Dict{Symbol,Float64}(), String[]),
                    false,
                    :snapshot,
                ),
            )
        end
        Uᴴ = Unew
        nprevious = n̄
        previous_scba = scba
        Jprevious = J
        populations_previous = populations
        if convergence_streak >= options.convergence.required_consecutive_poisson_passes
            outer_converged = true
            status = :converged
            break
        end
    end
    # No implicit extra SCBA solve: the latest solved snapshot is already
    # aligned. Explicit certify_solution owns any additional polish budget.
    previous_scba === nothing &&
        throw(ArgumentError("no SCBA state available for final audit"))
    final_scba = previous_scba
    Uᴴ = published_U
    nfinal = _electron_density_bar(problem, final_scba.green.Gˡ)
    observables = _collect_observables(problem, final_scba, Uᴴ)
    observables[:warnings] = warnings
    observables[:acceptance_measurement_context] = Dict{String,Any}(
        "backend"=>"production",
        "hilbert_roundoff_measured"=>production_options.verify_fft_roundoff || all(
            ==(:electron_electron),
            problem.kernels.enabled,
        ),
        "fresh_candidate_required"=>true,
    )
    observables[:inner_working_policy] = Dict{String,Any}(
        "mode"=>String(options.convergence.mode),
        "evaluations_since_resume"=>inner_policy_history,
        "strict_final_required"=>true,
        "forcing_rule"=>"min declared coarse band and 0.1 * best finite max(raw r_U, raw r_n); strict within 10 outer tolerances",
    )
    observables[:restart_contract] = restart_contract
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
    candidate = nothing
    audit_failure = final_scba.quality === :invalid ? String(final_scba.status) : nothing
    if audit_failure === nothing
        audit_field_bytes =
            sizeof(ComplexF64)*problem.numerical.N_E*problem.numerical.N_k*problem.numerical.N_b^2
        audit_phase = _production_phase_begin(
            production_options,
            :final_audit,
            0;
            stage = :audit,
            task_width = _production_thread_worker_count(
                production_options,
                problem.numerical.N_E*problem.numerical.N_k,
            ),
            memory_burst_bytes = (18+3length(problem.kernels.enabled))*audit_field_bytes,
        )
        audit_status = :failed
        try
            audit_sc, audit_roundoff = _scattering_candidate_production(
                problem,
                final_scba.green,
                cache,
                audit_phase.options;
                return_roundoff = true,
            )
            audit_embedding, audit_plus, audit_minus = _embedding_self_energy_selected(
                problem,
                final_scba.green,
                cache,
                audit_phase.options;
                scattering = audit_sc,
                h = project_hamiltonians(problem, Uᴴ),
            )
            candidate = FixedPointAuditCandidate(
                audit_sc,
                audit_embedding,
                audit_plus,
                audit_minus,
                audit_roundoff,
            )
            audit_status = :completed
        catch error
            error isa DomainError || rethrow()
            audit_failure = sprint(showerror, error)
        finally
            _production_phase_end(production_options, audit_phase; status = audit_status)
        end
    end
    quality_report =
        _final_quality_report(provisional, outer_converged; candidate, audit_failure)
    _publish_charge_diagnostics!(observables, quality_report.report)
    _publish_model_diagnostics!(observables, problem, final_scba, Uᴴ)
    final_warning = _final_validation_warning(provisional, quality_report)
    final_warning === nothing || push!(warnings, final_warning)
    report = quality_report.report
    converged = quality_report.converged
    final_status = if converged
        :converged
    elseif quality_report.approximate
        :approximate
    elseif final_scba.status === :research_continue && !outer_converged
        :outer_limit_with_warning
    elseif !final_scba.converged && outer_converged
        :final_scba_failed
    elseif !final_scba.converged && !startswith(String(status), "scba_")
        Symbol(status, :_final_scba_, final_scba.status)
    elseif outer_converged && !report.passed
        :validation_failed
    else
        status
    end
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
    state_observer === nothing || state_observer(result)
    checkpoint_sink === nothing || checkpoint_sink(result)
    _end_solver_stage(
        production_options,
        :poisson,
        owns_poisson_stage;
        status = (converged || quality_report.approximate) ? :completed : :incomplete,
        iteration = length(outer_history),
        total = options.max_poisson,
        metrics = SolverMetric(
            :J,
            _solver_current_density(
                Float64(ustrip(u"A/m^2", result.observables[:electron_flow_current])),
            ),
            "A/cm^2",
        ),
        message = String(final_status),
    )
    return result
end
