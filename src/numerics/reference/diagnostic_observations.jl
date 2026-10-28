"""
    solution_metric_diagnostics(solution)

Read detached summary values and strict/optional approximate reporting limits
from the stored state. Final report values override history values. No fresh map
is evaluated. Use `solution_scientific_assessment` for measurement and acceptance
status; the presence of a limit alone never establishes a measured pass.
"""
function solution_metric_diagnostics(solution::NEGFSolution)
    metrics = solution.report.metrics
    approximate = if solution.options.convergence.diagnostic_quality.enabled
        _solution_validation_limits(_approximate_options(solution.options).tolerances, metrics)
    else
        nothing
    end
    return SolutionMetricDiagnostics(
        _solution_report_metrics(solution),
        _solution_validation_limits(solution.options.tolerances, metrics),
        approximate,
    )
end

const _DIAGNOSTIC_MATRIX_KINDS =
    (:spectral, :occupied, :unoccupied, :broadening, :raw_occupied, :raw_unoccupied)

function _public_positivity_witness(w::_PositivityWitness)
    measured = w.matrix_kind !== :none && w.energy_index > 0 && w.momentum_index > 0 &&
        all(isfinite, (w.minimum_eigenvalue, w.block_norm, w.floor, w.ratio,
            w.absolute_defect, w.relative_defect, w.backward_error,
            w.hermiticity_defect, w.construction_scale))
    return PositivityWitness(
        measured ? :measured : :error,
        (getfield(w, name) for name in fieldnames(_PositivityWitness))...,
    )
end

function _diagnostic_matrix_eigenvalues(block)
    all(isfinite, block) || return fill(NaN, size(block, 1))
    return try
        eigvals(Hermitian((block + block') / 2))
    catch error
        error isa LinearAlgebra.LAPACKException || rethrow()
        fill(NaN, size(block, 1))
    end
end

function _check_diagnostic_state_shapes(solution::NEGFSolution)
    problem = solution.problem
    n = problem.numerical
    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    n.N_E > 0 && n.N_k > 0 && n.N_b > 0 ||
        throw(DimensionMismatch("diagnostic axes must be nonempty"))
    for field in (:Gᴿ, :Gˡ, :Gᵍ, :A)
        size(getfield(solution.scba.green, field)) == shape ||
            throw(DimensionMismatch("diagnostic Green shape differs from problem"))
    end
    for family in (solution.scba.embedding,
        (solution.scba.scattering[name] for name in problem.kernels.enabled)...)
        for field in (:Σᴿ, :Σˡ, :Σᵍ)
            size(getfield(family, field)) == shape ||
                throw(DimensionMismatch("diagnostic self-energy shape differs from problem"))
        end
    end
    length(problem.grids.wᴱ) == n.N_E && length(problem.grids.wᵏ) == n.N_k ||
        throw(DimensionMismatch("diagnostic quadrature axes differ from problem"))
    return nothing
end

# Numerical body formerly owned by Runner's final matrix audit. Its k→E scan,
# quadrature order, strict tie comparison and failed IEEE payload are retained.
function _matrix_integral_diagnostics(solution::NEGFSolution)
    problem, green = solution.problem, solution.scba.green
    NE, Nk, Nb = problem.numerical.N_E, problem.numerical.N_k, problem.numerical.N_b
    integral = zeros(ComplexF64, Nk, Nb, Nb)
    cumulative = zeros(Float64, NE, Nk)
    eigenvalues = zeros(Float64, Nk, Nb)
    negative = zeros(Float64, 6, Nk)
    negative_blocks, beyond_budget, invalid_blocks = zeros(Int, 6), zeros(Int, 6), zeros(Int, 6)
    witnesses = fill(_empty_positivity_witness(), 6)
    for k = 1:Nk
        for e = 1:NE
            weight = problem.grids.wᴱ[e] / (2π)
            integral[k, :, :] .+= weight .* view(green.A, e, k, :, :)
            cumulative[e, k] = real(tr(view(integral, k, :, :)))
            blocks = _scba_positivity_blocks(solution.scba, problem.kernels.enabled, e, k)
            for (index, (kind, block)) in enumerate(blocks)
                values = _diagnostic_matrix_eigenvalues(block)
                witness = _positivity_witness(kind, e, k, values;
                    hermiticity_defect = all(isfinite, block) ? norm(block - block') / 2 : Inf)
                witness.ratio > witnesses[index].ratio && (witnesses[index] = witness)
                if !all(isfinite, values)
                    negative[index, k] = NaN
                    invalid_blocks[index] += 1
                    continue
                end
                negative[index, k] += weight * sum(value -> max(0.0, -value), values)
                negative_blocks[index] += minimum(values) < 0
                beyond_budget[index] += minimum(values) < -psd_error_budget(maximum(abs, values), Nb)
            end
        end
        eigenvalues[k, :] = _diagnostic_matrix_eigenvalues(Matrix(view(integral, k, :, :)))
    end
    phase_space_weights = ntuple(5) do j
        i = (1, 2, 3, 5, 6)[j]
        problem.physical.g_s * sum(negative[i, :] .* problem.grids.wᵏ) / problem.scales.L₀_m^2
    end
    integrals = MatrixIntegralDiagnostics(integral, eigenvalues, cumulative, negative,
        Tuple(negative_blocks), Tuple(beyond_budget), Tuple(invalid_blocks), NE * Nk,
        phase_space_weights)
    return witnesses, integrals
end

"""
    positivity_diagnostics(solution; integrate=false)

Observe all six PSD matrix kinds without clipping or changing the state. The
light E→k scan selects worst witnesses. `integrate=true` instead performs the
historical k→E matrix-audit scan and returns full integrals/counts/negative
weights, reusing its witnesses. Both scans keep first-on-tie (`>`) selection.
Local `< r_PSD` flags do not establish global scientific acceptance. Invalid
blocks retain IEEE fields and `status=:error`; omitted integrals are `nothing`.
"""
function positivity_diagnostics(solution::NEGFSolution; integrate::Bool = false)
    _check_diagnostic_state_shapes(solution)
    private_witnesses, integrals = if integrate
        _matrix_integral_diagnostics(solution)
    else
        (_scba_positivity_witnesses(solution.scba, solution.problem.kernels.enabled), nothing)
    end
    witnesses = ntuple(i -> _public_positivity_witness(private_witnesses[i]), 6)
    limit = solution.options.tolerances.r_PSD
    absolute = ntuple(i -> witnesses[i].backward_error + limit * witnesses[i].block_norm, 6)
    satisfied = ntuple(i -> witnesses[i].status === :measured && witnesses[i].ratio < limit, 6)
    return PositivityDiagnostics(PSD_METRIC_VERSION, _DIAGNOSTIC_MATRIX_KINDS, witnesses,
        limit, absolute, satisfied, integrate ? :momentum_energy : :energy_momentum, integrals)
end

"""
    basis_validation_diagnostics(basis, projected, numerical)

Measure the fixed-basis subset of `validate_problem` using the same limits:
orthogonality/eigen residuals, T-adjoint, projected Hermiticity and conditional
PzP translation metrics. `projected` has `(Nk,Nb,Nb)` axes. No scattering problem
or transport state is built; this subset is not full static validation.
"""
function basis_validation_diagnostics(basis::BasisData,
    projected::AbstractArray{<:Number,3}, n::NumericalParameters)
    size(projected) == (n.N_k, n.N_b, n.N_b) && n.N_k > 0 && n.N_b > 0 ||
        throw(DimensionMismatch("projected Hamiltonian must have (Nk,Nb,Nb) axes"))
    size(basis.T₊) == size(basis.T₋) == (n.N_b, n.N_b) ||
        throw(DimensionMismatch("basis translations must have (Nb,Nb) axes"))
    metrics = Dict(:basis_orthogonality => basis.r_orth, :basis_eigen => basis.r_eigen,
        :T_adjoint => norm(basis.T₋ - basis.T₊') / (norm(basis.T₊) + 1e-14),
        :hamiltonian_hermiticity => maximum(
            norm(Matrix(view(projected, m, :, :)) - Matrix(view(projected, m, :, :))') /
            (norm(Matrix(view(projected, m, :, :))) + 1e-14) for m = 1:n.N_k))
    limits = Dict(name => _STRUCTURAL_VALIDATION_LIMITS[name] for name in keys(metrics))
    if basis.localization === :pzp
        metrics[:translation] = basis.r_translation
        metrics[:translation_nearest] = basis.r_translation_nearest
        metrics[:translation_two_pairs] = basis.r_translation_two_pairs
        limits[:translation] = _PZP_TRANSLATION_LIMIT
        limits[:translation_nearest] = _PZP_TRANSLATION_LIMIT
        limits[:translation_two_pairs] = n.P_basis >= 3 ? _PZP_TWO_PAIR_LIMIT : _PZP_TRANSLATION_LIMIT
    end
    messages = String[]
    for (name, limit) in limits
        (!isfinite(metrics[name]) || metrics[name] > limit) &&
            push!(messages, "$name=$(metrics[name]) exceeds $limit")
    end
    return BasisValidationDiagnostics(ConvergenceReport(isempty(messages), metrics, messages), limits)
end

"""
    projected_hamiltonian_spectrum(h)

Evaluate one diagnostic projected eigensystem using the existing effective-level
Hermiticity gate and floor. Reject nonsquare, nonfinite or non-Hermitian input;
do not repair the matrix. Input energy units are retained.
"""
function projected_hamiltonian_spectrum(h::AbstractMatrix{<:Number})
    size(h, 1) == size(h, 2) && size(h, 1) > 0 ||
        throw(DimensionMismatch("projected Hamiltonian must be nonempty and square"))
    all(isfinite, h) || throw(DomainError(h, "projected Hamiltonian is nonfinite"))
    floor = _EFFECTIVE_HAMILTONIAN_NORM_FLOOR
    limit = _EFFECTIVE_HAMILTONIAN_HERMITICITY_LIMIT
    residual = norm(h - h') / (norm(h) + floor)
    residual <= limit || throw(DomainError(residual, "projected h(k=0) is not Hermitian"))
    decomposition = eigen(Hermitian(h))
    values, vectors = Vector{Float64}(decomposition.values), Matrix{ComplexF64}(decomposition.vectors)
    eigen_residual = norm(h * vectors - vectors * Diagonal(values)) / (norm(h) + floor)
    return ProjectedHamiltonianSpectrum(values, vectors, residual, eigen_residual, floor, limit)
end
