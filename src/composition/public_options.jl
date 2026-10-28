"""
    validate_options(options)

Validate one SolverOptions, AlgorithmOptions, ProductionOptions or
ProductionKernelOptions value using its owning admission rules. Return the same
value or throw the existing ArgumentError; do not build a problem or call sinks.
This establishes option coherence, not scientific validity.
"""
validate_options(options::SolverOptions) = _check_options(options)
validate_options(options::AlgorithmOptions) = _check_algorithm_options(options)
validate_options(options::ProductionOptions) = _check_production_options(options)
validate_options(options::ProductionKernelOptions) = _check_production_kernel_options(options)

"""
    validate_result_classification(converged, status, scba_quality)

Check coherence of the three result labels and return `scba_quality`. Strict
inner SCBA quality may coexist with unconverged outer iteration. This label
check does not measure or attest scientific acceptance.
"""
validate_result_classification(converged::Bool, status::Symbol, scba_quality::Symbol) =
    _check_scba_result_classification(converged, status, scba_quality)
