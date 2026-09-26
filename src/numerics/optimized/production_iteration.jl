"""Measure the fresh raw fixed-point map before mixing, with reusable residual storage."""
function _production_fixed_point_measurement(
    problem::NEGFProblem,
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
    lambda,
    workspace;
    worker_count::Integer = 0,
)
    residuals = production_residual_suite(
        problem.grids.ε,
        h,
        green,
        total,
        candidate_total;
        workspace = workspace,
        worker_count,
    )
    rD = residuals.r_D
    rA = residuals.r_A
    rK = residuals.r_K
    rPSD = residuals.r_PSD
    rcaus = residuals.r_caus
    rΣ = max(
        production_selfenergy_residual(
            embedding,
            candidate_embedding;
            workspace = workspace,
            worker_count,
        ),
        production_selfenergy_residual(
            plus,
            candidate_plus;
            workspace = workspace,
            worker_count,
        ),
        production_selfenergy_residual(
            minus,
            candidate_minus;
            workspace = workspace,
            worker_count,
        ),
    )
    for mechanism in problem.kernels.enabled
        rΣ = max(
            rΣ,
            production_selfenergy_residual(
                scattering[mechanism],
                candidate_sc[mechanism];
                workspace = workspace,
                worker_count,
            ),
        )
    end
    rλ = abs(lambda - 1)
    return (; rD, rA, rK, rΣ, rλ, rPSD, rcaus)
end

"""Exact block trace in O(N_b²), without materializing either matrix product.

The educational observable remains the independent matrix-product oracle.
Each diagonal dot product is accumulated separately before the trace sum.
"""
function _production_current_trace(Σᵍ, Σˡ, Gˡ, Gᵍ, e, m)
    Nb = size(Gˡ, 3)
    value = 0.0 + 0.0im
    @inbounds for a = 1:Nb
        outgoing = 0.0 + 0.0im
        incoming = 0.0 + 0.0im
        for b = 1:Nb
            outgoing += Σᵍ[e, m, a, b] * Gˡ[e, m, b, a]
            incoming += Σˡ[e, m, a, b] * Gᵍ[e, m, b, a]
        end
        value += outgoing - incoming
    end
    return value
end

function _production_boundary_flux_bar(problem, family, green)
    # Preserve the reference cancellation-aware energy accumulation order.
    energy = problem.grids.ε
    order = sortperm(collect(eachindex(energy)); by = e -> (abs(energy[e]), energy[e]))
    value = 0.0 + 0.0im
    for e in order, m in axes(green.Gᴿ, 2)
        value +=
            problem.grids.wᴱ[e] *
            problem.grids.wᵏ[m] *
            _production_current_trace(family.Σᵍ, family.Σˡ, green.Gˡ, green.Gᵍ, e, m)
    end
    return real(problem.physical.g_s * value / (2π))
end

"""Accumulate only the diagonal populations needed by the SCBA change gate."""
function _production_state_populations(problem, Gˡ)
    populations = zeros(ComplexF64, size(Gˡ, 3))
    prefactor = -im * problem.physical.g_s / (2π)
    @inbounds for e in axes(Gˡ, 1), m in axes(Gˡ, 2)
        weight = prefactor * problem.grids.wᴱ[e] * problem.grids.wᵏ[m]
        for a in eachindex(populations)
            populations[a] += weight * Gˡ[e, m, a, a]
        end
    end
    target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    return real.(populations) ./ target
end

"""Scientific observable changes of the current Green state, independent of mixing strength."""
function _production_observable_measurement(
    problem::NEGFProblem,
    green,
    plus,
    Jprevious,
    populations_previous,
)
    Jbar = _production_boundary_flux_bar(problem, plus, green)
    J = Jbar * Float64(ustrip(u"A/m^2", problem.scales.J₀))
    target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    populations = _production_state_populations(problem, green.Gˡ)
    rJchange =
        Jprevious === nothing ? Inf :
        abs(J - Jprevious) / max(abs(J), abs(Jprevious), 1e-12)
    rpopulation =
        populations_previous === nothing ? Inf :
        norm(populations - populations_previous) / (norm(populations) + 1e-14)
    return (; J, target, populations, rJchange, rpopulation)
end

"""Reconstruct one measured PSD witness from this iteration, never a later state.

All blocks use dimensionless solver units. Raw direct Gp and the reconstructed
A-Gn are deliberately both retained to separate input, subtraction and number-
constraint errors. Per-channel in/out blocks support critical-vector projections.
"""
function _scba_physics_witness(witness, green, total, scattering, plus, minus)
    e, m = witness.energy_index, witness.momentum_index
    e > 0 && m > 0 || return nothing
    GR = _matrix_block(green.Gᴿ, e, m)
    sigmaL = _matrix_block(total.Σˡ, e, m)
    sigmaG = _matrix_block(total.Σᵍ, e, m)
    A = _matrix_block(green.A, e, m)
    Gn = -im .* _matrix_block(green.Gˡ, e, m)
    Gp = im .* _matrix_block(green.Gᵍ, e, m)
    matrices = Dict{String,Matrix{ComplexF64}}(
        "GR_dimensionless" => GR,
        "A_dimensionless" => A,
        "Gn_raw_dimensionless" => -im .* (GR*sigmaL*GR'),
        "Gp_raw_direct_dimensionless" => im .* (GR*sigmaG*GR'),
        "Gn_normalized_dimensionless" => Gn,
        "Gp_normalized_dimensionless" => Gp,
        "Gp_reconstructed_dimensionless" => A-Gn,
        "SigmaR_total_dimensionless" => _matrix_block(total.Σᴿ, e, m),
        "Sigma_in_total_dimensionless" => -im .* sigmaL,
        "Sigma_out_total_dimensionless" => im .* sigmaG,
        "Gamma_dimensionless" => im .* (sigmaG-sigmaL),
    )
    for (name, family) in ((String(k), v) for (k, v) in scattering)
        matrices["Sigma_in_$(name)_dimensionless"] = -im .* _matrix_block(family.Σˡ, e, m)
        matrices["Sigma_out_$(name)_dimensionless"] = im .* _matrix_block(family.Σᵍ, e, m)
    end
    for (name, family) in (("embedding_plus", plus), ("embedding_minus", minus))
        matrices["Sigma_in_$(name)_dimensionless"] = -im .* _matrix_block(family.Σˡ, e, m)
        matrices["Sigma_out_$(name)_dimensionless"] = im .* _matrix_block(family.Σᵍ, e, m)
    end
    name =
        witness.matrix_kind === :spectral ? "A_dimensionless" :
        witness.matrix_kind === :occupied ? "Gn_normalized_dimensionless" :
        witness.matrix_kind === :unoccupied ? "Gp_normalized_dimensionless" :
        witness.matrix_kind === :raw_occupied ? "Gn_raw_dimensionless" :
        witness.matrix_kind === :raw_unoccupied ? "Gp_raw_direct_dimensionless" :
        "Gamma_dimensionless"
    critical = matrices[name]
    vector =
        all(isfinite, critical) ?
        Vector{ComplexF64}(eigen(Hermitian((critical+critical')/2)).vectors[:, 1]) :
        ComplexF64[]
    return SCBAPhysicsWitness(
        witness.matrix_kind,
        e,
        m,
        witness.minimum_eigenvalue,
        witness.block_norm,
        witness.backward_error,
        witness.absolute_defect,
        witness.relative_defect,
        witness.ratio,
        witness.hermiticity_defect,
        matrices,
        vector,
    )
end

"""Bound in-memory matrix payloads while retaining every scalar iteration witness."""
mutable struct _SCBAWitnessRetention
    first::Int
    worst::Int
    last::Int
end

function _SCBAWitnessRetention(history)
    indices = findall(row -> row.witness !== nothing, history)
    isempty(indices) && return _SCBAWitnessRetention(0, 0, 0)
    worst = first(indices)
    for index in indices
        history[index].witness.ratio > history[worst].witness.ratio && (worst=index)
    end
    return _SCBAWitnessRetention(first(indices), worst, last(indices))
end

function _strip_scba_witness_payload(row::SCBAIteration)
    w = row.witness
    w === nothing && return row
    isempty(w.matrices) && isempty(w.eigenvector) && return row
    scalar = SCBAPhysicsWitness(
        w.matrix_kind,
        w.energy_index,
        w.momentum_index,
        w.minimum_eigenvalue,
        w.block_norm,
        w.backward_error,
        w.absolute_defect,
        w.relative_defect,
        w.ratio,
        w.hermiticity_defect,
        Dict{String,Matrix{ComplexF64}}(),
        ComplexF64[],
    )
    return SCBAIteration(
        (
            name === :witness ? scalar : getfield(row, name) for
            name in fieldnames(SCBAIteration)
        )...,
    )
end

function _retain_scba_witness_payloads!(history, retention::_SCBAWitnessRetention)
    index = length(history)
    history[index].witness === nothing && return nothing
    previous_worst, previous_last = retention.worst, retention.last
    retention.first == 0 && (retention.first=index)
    if retention.worst == 0 ||
       history[index].witness.ratio > history[retention.worst].witness.ratio
        retention.worst=index
    end
    retention.last=index
    for old in (previous_worst, previous_last)
        old == 0 && continue
        old in (retention.first, retention.worst, retention.last) && continue
        history[old] = _strip_scba_witness_payload(history[old])
    end
    return nothing
end
