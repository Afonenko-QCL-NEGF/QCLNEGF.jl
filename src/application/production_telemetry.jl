function _emit_scba_progress(
    options::ProductionOptions,
    row::SCBAIteration,
    total::Int,
    assessments,
    timings,
    quality::Symbol,
)
    iteration = row.ν
    (; assessment, convergence_streak, diagnostic, diagnostic_quality_streak,
       fixed_point, positivity_witness) = assessments
    if options.event_sink !== nothing
        metrics = SolverMetric[
            SolverMetric(:J, _solver_current_density(row.J), "A/cm^2"),
            SolverMetric(:r_D, row.r_D),
            SolverMetric(:r_A, row.r_A),
            SolverMetric(:r_K, row.r_K),
            SolverMetric(:r_Σ, row.r_Σ),
            SolverMetric(:r_λ, row.r_λ),
            SolverMetric(:r_PSD, row.r_PSD),
            SolverMetric(:fixed_point_limiting_gate, String(fixed_point.limiting_metric)),
            SolverMetric(:fixed_point_limiting_ratio, fixed_point.limiting_ratio),
            SolverMetric(:fixed_point_passed, fixed_point.passed),
            SolverMetric(:psd_matrix_kind, String(positivity_witness.matrix_kind)),
            SolverMetric(:psd_energy_index, positivity_witness.energy_index),
            SolverMetric(:psd_momentum_index, positivity_witness.momentum_index),
            SolverMetric(:psd_minimum_eigenvalue_scaled, positivity_witness.minimum_eigenvalue),
            SolverMetric(:psd_block_norm_scaled, positivity_witness.block_norm),
            SolverMetric(:psd_floor_scaled, positivity_witness.floor),
            SolverMetric(:psd_metric_version, PSD_METRIC_VERSION),
            SolverMetric(:psd_absolute_defect, positivity_witness.absolute_defect),
            SolverMetric(:psd_relative_defect, positivity_witness.relative_defect),
            SolverMetric(:psd_backward_error, positivity_witness.backward_error),
            SolverMetric(:psd_hermiticity_defect, positivity_witness.hermiticity_defect),
            SolverMetric(:r_caus, row.r_caus),
            SolverMetric(:r_roundoff, row.r_roundoff),
            SolverMetric(:r_Jchange, row.r_Jchange),
            SolverMetric(:r_population, row.r_population),
            SolverMetric(:limiting_gate, String(assessment.limiting_metric)),
            SolverMetric(:limiting_ratio, assessment.limiting_ratio),
            SolverMetric(:failed_gate_count, length(assessment.failed_metrics)),
            SolverMetric(:convergence_streak, convergence_streak),
            SolverMetric(:diagnostic_quality, String(quality)),
            SolverMetric(:outer_iteration, options.outer_iteration),
            SolverMetric(
                :diagnostic_quality_enabled,
                diagnostic.eligible ||
                    diagnostic.limiting_metric !== :diagnostic_quality_disabled,
            ),
            SolverMetric(:diagnostic_quality_streak, diagnostic_quality_streak),
            SolverMetric(:diagnostic_limiting_gate, String(diagnostic.limiting_metric)),
            SolverMetric(:diagnostic_limiting_ratio, diagnostic.limiting_ratio),
            SolverMetric(:t_total, timings.total, "s"),
        ]
        if options.phase_timing
            append!(
                metrics,
                SolverMetric[
                    SolverMetric(:t_candidate, timings.candidate, "s"),
                    SolverMetric(:t_candidate_scattering, timings.candidate_scattering, "s"),
                    SolverMetric(:t_candidate_embedding, timings.candidate_embedding, "s"),
                    SolverMetric(:t_candidate_assembly, timings.candidate_assembly, "s"),
                    SolverMetric(:t_physics_markers, timings.physics_markers, "s"),
                    SolverMetric(:t_dyson, timings.dyson, "s"),
                    SolverMetric(:t_mixing, timings.mixing, "s"),
                    SolverMetric(:t_observables, timings.observables, "s"),
                    SolverMetric(:t_residuals, timings.residual, "s"),
                ],
            )
        end
        _update_solver_stage(options, :scba; iteration, total, metrics) && return nothing
    end
    every = options.progress_every_scba
    if every > 0 && (iteration == 1 || iteration % every == 0 || iteration == total)
        @info "SCBA progress" iteration total rD=row.r_D rK=row.r_K rΣ=row.r_Σ rλ=row.r_λ current_A_per_cm2=_solver_current_density(row.J) total_seconds=timings.total
    end
    return nothing
end

"""Снимок времени и счётчиков CPU, сборки мусора и выделений по всему процессу."""
function _production_counter_sample()
    gc = Base.gc_num()
    return (
        monotonic_ns = time_ns(),
        cpu_seconds = Float64(ccall(:clock, Clong, ())) / 1.0e6,
        allocated_bytes = Int64(Base.gc_bytes()),
        gc_seconds = Float64(gc.total_time) * 1e-9,
        gc_count = Int64(gc.pause),
    )
end

function _production_phase_begin(
    options::ProductionOptions,
    name::Symbol,
    iteration::Int;
    task_width::Int = 1,
    work_units::Int = 0,
    stage::Symbol = :scba,
    elastic::Bool = true,
    workspace_bytes::Union{Nothing,Int} = nothing,
    memory_burst_bytes::Union{Nothing,Int} = nothing,
)
    if options.phase_request !== nothing
        granted = options.phase_request((;
            name,
            stage,
            iteration,
            task_width,
            work_units,
            elastic,
            workspace_bytes,
            memory_burst_bytes,
            safe_boundary = true,
        ))
        granted isa Integer && 1 <= granted <= task_width ||
            throw(ArgumentError("phase task width exceeds the declared task allocation"))
        !elastic &&
            granted != task_width &&
            throw(ArgumentError("fixed-workspace phase cannot change its worker width"))
        task_width = Int(granted)
    end
    phase_options = with_production_options(options; worker_count = task_width)
    recording = options.phase_timing && options.event_sink !== nothing
    sample = recording ? _production_counter_sample() : nothing
    recording || return (;
        name,
        iteration,
        stage,
        task_width,
        work_units,
        sample,
        recording,
        workspace_bytes,
        memory_burst_bytes,
        options = phase_options,
    )
    metrics = SolverMetric[
        SolverMetric(:producer_monotonic_ns, sample.monotonic_ns, "ns"),
        SolverMetric(:start_monotonic_ns, sample.monotonic_ns, "ns"),
        SolverMetric(:cpu_seconds_total, sample.cpu_seconds, "s"),
        SolverMetric(:allocated_bytes_total, sample.allocated_bytes, "bytes"),
        SolverMetric(:gc_seconds_total, sample.gc_seconds, "s"),
        SolverMetric(:gc_count_total, sample.gc_count),
        SolverMetric(:task_width, task_width),
        SolverMetric(:work_units, work_units),
        SolverMetric(:counter_scope, "process"),
        SolverMetric(:outer_iteration, options.outer_iteration),
    ]
    push!(
        metrics,
        SolverMetric(
            :memory_estimate_status,
            memory_burst_bytes === nothing ? "unavailable" : "dimension_estimate",
        ),
    )
    workspace_bytes === nothing ||
        push!(metrics, SolverMetric(:workspace_bytes, workspace_bytes, "bytes"))
    memory_burst_bytes === nothing ||
        push!(metrics, SolverMetric(:memory_burst_bytes, memory_burst_bytes, "bytes"))
    options.event_sink(
        SolverEvent(
            :phase_begin,
            stage,
            String(name),
            :running,
            iteration,
            nothing,
            metrics,
            "",
        ),
    )
    return (;
        name,
        iteration,
        stage,
        task_width,
        work_units,
        sample,
        recording,
        workspace_bytes,
        memory_burst_bytes,
        options = phase_options,
    )
end

function _production_phase_end(
    options::ProductionOptions,
    token;
    status::Symbol = :completed,
)
    (token === nothing || !token.recording) && return nothing
    sample = _production_counter_sample()
    metrics = SolverMetric[
        SolverMetric(:producer_monotonic_ns, sample.monotonic_ns, "ns"),
        SolverMetric(:start_monotonic_ns, token.sample.monotonic_ns, "ns"),
        SolverMetric(:end_monotonic_ns, sample.monotonic_ns, "ns"),
        SolverMetric(
            :duration_seconds,
            (sample.monotonic_ns-token.sample.monotonic_ns)*1e-9,
            "s",
        ),
        SolverMetric(:cpu_seconds_total, sample.cpu_seconds, "s"),
        SolverMetric(:allocated_bytes_total, sample.allocated_bytes, "bytes"),
        SolverMetric(:gc_seconds_total, sample.gc_seconds, "s"),
        SolverMetric(:gc_count_total, sample.gc_count),
        SolverMetric(
            :cpu_seconds,
            max(0.0, sample.cpu_seconds-token.sample.cpu_seconds),
            "s",
        ),
        SolverMetric(
            :allocated_bytes,
            max(0, sample.allocated_bytes-token.sample.allocated_bytes),
            "bytes",
        ),
        SolverMetric(:gc_seconds, max(0.0, sample.gc_seconds-token.sample.gc_seconds), "s"),
        SolverMetric(:gc_count, max(0, sample.gc_count-token.sample.gc_count)),
        SolverMetric(:task_width, token.task_width),
        SolverMetric(:work_units, token.work_units),
        SolverMetric(:counter_scope, "process"),
        SolverMetric(:outer_iteration, options.outer_iteration),
    ]
    push!(
        metrics,
        SolverMetric(
            :memory_estimate_status,
            token.memory_burst_bytes === nothing ? "unavailable" : "dimension_estimate",
        ),
    )
    token.workspace_bytes === nothing ||
        push!(metrics, SolverMetric(:workspace_bytes, token.workspace_bytes, "bytes"))
    token.memory_burst_bytes === nothing ||
        push!(metrics, SolverMetric(:memory_burst_bytes, token.memory_burst_bytes, "bytes"))
    options.event_sink(
        SolverEvent(
            :phase_end,
            token.stage,
            String(token.name),
            status,
            token.iteration,
            nothing,
            metrics,
            "",
        ),
    )
    return nothing
end
