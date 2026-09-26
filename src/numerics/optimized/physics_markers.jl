# All routines in this file observe the current state. They never repair Green
# matrices, alter a self-energy, or participate in acceptance.
const _MARKER_ROUNDOFF_MULTIPLIER = 64
const _MARKER_FERMI_BRACKET_THERMAL_WIDTHS = 50
const _MARKER_LINEWIDTH_SAMPLING_METHOD = :spectral_midpoint_quantiles_A_eigenmodes_v1

function _marker_equilibrium_status(problem::NEGFProblem)
    iszero(problem.physical.F_bias) || return :finite_bias
    if problem.scattering.LO
        problem.models.lo_population === :thermal || return :nonthermal_LO_bath
        problem.physical.Tᴸ == problem.physical.Tᴸᴼ || return :unequal_bath_temperatures
    end
    ee = problem.models.electron_electron
    ee.mode === :none ||
        ee.electron_temperature == problem.physical.Tᴸ ||
        return :unequal_bath_temperatures
    return :available
end

"""Fit one equilibrium chemical potential using the spectral number, not normalized Gn."""
function _marker_spectral_chemical_potential(energy, spectral_number_weight, target, kBT)
    all(isfinite, spectral_number_weight) ||
        return (; mu = NaN, status = :nonfinite_spectrum)
    floor = _MARKER_ROUNDOFF_MULTIPLIER*eps(Float64)*maximum(abs, spectral_number_weight)
    any(x -> x < -floor, spectral_number_weight) &&
        return (; mu = NaN, status = :negative_spectral_weight)
    # Ignore only floating-point sign noise in this diagnostic number fit.
    weights = max.(spectral_number_weight, 0.0)
    charge(mu) = sum(weights[e]*_fermi(energy[e], mu, kBT) for e in eachindex(energy))
    lo = first(energy) - _MARKER_FERMI_BRACKET_THERMAL_WIDTHS*kBT
    hi = last(energy) + _MARKER_FERMI_BRACKET_THERMAL_WIDTHS*kBT
    charge(lo) <= target <= charge(hi) ||
        return (; mu = NaN, status = :spectral_number_not_bracketed)
    mu = find_zero(x -> charge(x)-target, (lo, hi), Bisection())
    return (; mu, status = :available)
end

function _marker_weighted_quantile(values, weights, probability)
    order = sortperm(values)
    threshold = probability*sum(weights)
    accumulated = 0.0
    for index in order
        accumulated += weights[index]
        accumulated >= threshold && return values[index]
    end
    return values[last(order)]
end

"""Spectrally weighted *sampled* modal Γ/ΔE, bounded by the policy block budget.

Select midpoint quantiles of wE*wK*Tr(A), then resolve each selected block in
its spectral eigenvectors. Its modal linewidth is u'Γu, weighted by that
mode's positive spectral weight. Γ has energy units, so no hbar or 2pi is
inserted. Overlapping resonances need not have this projected width as FWHM.
"""
function _sampled_linewidth_markers(A, total, grids, policy::SCBAPhysicsMarkerPolicy)
    NE, Nk, Nb, _ = size(A)
    step = (last(grids.ε)-first(grids.ε))/(NE-1)
    missing(status) =
        (; q10 = NaN, q50 = NaN, q90 = NaN, underresolved = NaN, blocks = 0, status)
    isfinite(step) && step > 0 || return missing(:invalid_energy_grid)
    weights = zeros(Float64, NE*Nk)
    @inbounds for m = 1:Nk, e = 1:NE
        weights[e+(m-1)*NE] = grids.wᴱ[e]*grids.wᵏ[m]*sum(real(A[e, m, a, a]) for a = 1:Nb)
    end
    all(isfinite, weights) || return missing(:nonfinite_spectrum)
    roundoff = _MARKER_ROUNDOFF_MULTIPLIER*Nb*eps(Float64)*maximum(abs, weights)
    any(x -> x < -roundoff, weights) && return missing(:negative_spectral_weight)
    weights .= max.(weights, 0.0)
    mass = sum(weights)
    mass > 0 || return missing(:empty_spectrum)
    cumulative = cumsum(weights)
    count = min(policy.max_spectral_blocks, length(weights))
    selected = [
        min(searchsortedfirst(cumulative, (j-0.5)*mass/count), length(weights)) for
        j = 1:count
    ]
    values, mode_weights = Float64[], Float64[]
    blocks = 0
    status = :available
    previous = 0
    for linear in selected
        # Repeated quantiles count their probability mass, but use one eigensolve.
        linear == previous && continue
        previous = linear
        multiplicity =
            searchsortedlast(selected, linear)-searchsortedfirst(selected, linear)+1
        e, m = (linear-1)%NE+1, (linear-1)÷NE+1
        block = _matrix_block(A, e, m)
        all(isfinite, block) || return missing(:nonfinite_spectrum)
        spectral = eigen(Hermitian((block+block')/2))
        scale = maximum(abs, spectral.values)
        tolerance = _MARKER_ROUNDOFF_MULTIPLIER*Nb*eps(Float64)*scale
        minimum(spectral.values) < -tolerance && (status = :indefinite_spectrum)
        gamma = im .* (_matrix_block(total.Σᵍ, e, m)-_matrix_block(total.Σˡ, e, m))
        all(isfinite, gamma) || return missing(:nonfinite_broadening)
        gamma = (gamma+gamma')/2
        gamma_floor = _MARKER_ROUNDOFF_MULTIPLIER*Nb*eps(Float64)*norm(gamma)
        cutoff = policy.relative_mode_weight_floor*max(maximum(spectral.values), 0.0)
        kept = findall(x -> x > cutoff, spectral.values)
        isempty(kept) && continue
        normalization = sum(spectral.values[index] for index in kept)
        for index in kept
            vector = @view spectral.vectors[:, index]
            linewidth = real(dot(vector, gamma*vector))
            linewidth < -gamma_floor && (status = :noncausal_broadening)
            push!(values, linewidth/step)
            push!(mode_weights, multiplicity*spectral.values[index]/normalization)
        end
        blocks += 1
    end
    isempty(values) && return missing(:empty_spectrum)
    return (;
        q10 = _marker_weighted_quantile(values, mode_weights, 0.1),
        q50 = _marker_weighted_quantile(values, mode_weights, 0.5),
        q90 = _marker_weighted_quantile(values, mode_weights, 0.9),
        underresolved = sum(
            mode_weights[i] for i in eachindex(values) if values[i] < 1;
            init = 0.0,
        )/sum(mode_weights),
        blocks,
        status,
    )
end

_marker_relative_norm(error_squared, scale_squared) =
    scale_squared > 0 ? sqrt(error_squared/scale_squared) :
    iszero(error_squared) ? 0.0 : Inf

"""No array-sized temporaries: the diagnostic retains the unmixed map scale."""
function _channel_component_markers!(rows, channel, mixed, fresh)
    for (component, field) in ((:retarded, :Σᴿ), (:lesser, :Σˡ), (:greater, :Σᵍ))
        old, new = getfield(mixed, field), getfield(fresh, field)
        size(old) == size(new) || throw(DimensionMismatch("marker self-energy shape"))
        difference, scale = 0.0, 0.0
        @inbounds for index in eachindex(old, new)
            difference += abs2(new[index]-old[index])
            scale += abs2(new[index])
        end
        absolute, denominator = sqrt(difference), sqrt(scale)
        push!(
            rows,
            SCBAChannelMarker(
                channel,
                component,
                absolute,
                denominator,
                absolute/(denominator+1e-14),
            ),
        )
    end
    return rows
end

"""Compute particle/energy balances on exactly the same G without block copies."""
function _channel_collision_marker(problem, green, channel, state_kind, family)
    particle, energy = 0.0im, 0.0im
    particle_absolute, energy_absolute = 0.0, 0.0
    NE, Nk, Nb, _ = size(green.Gᴿ)
    @inbounds for e = 1:NE, m = 1:Nk
        incoming, outgoing = 0.0im, 0.0im
        for b = 1:Nb, a = 1:Nb
            incoming += family.Σˡ[e, m, a, b]*green.Gᵍ[e, m, b, a]
            outgoing += family.Σᵍ[e, m, a, b]*green.Gˡ[e, m, b, a]
        end
        weight = problem.grids.wᴱ[e]*problem.grids.wᵏ[m]
        epsilon = problem.grids.ε[e]
        particle += weight*(incoming-outgoing)
        energy += weight*epsilon*(incoming-outgoing)
        scale = weight*(abs(incoming)+abs(outgoing))
        particle_absolute += scale
        energy_absolute += abs(epsilon)*scale
    end
    prefactor = problem.physical.g_s/(2π)
    return SCBACollisionMarker(
        channel,
        state_kind,
        prefactor*real(particle),
        prefactor*particle_absolute,
        prefactor*imag(particle),
        prefactor*real(energy),
        prefactor*energy_absolute,
        prefactor*imag(energy),
    )
end

function _fresh_map_markers(
    problem,
    green,
    total;
    mixed_scattering = nothing,
    fresh_scattering = nothing,
    mixed_embedding = nothing,
    fresh_embedding = nothing,
    mixed_plus = nothing,
    fresh_plus = nothing,
    mixed_minus = nothing,
    fresh_minus = nothing,
    fresh_total = nothing,
)
    channels, collisions = SCBAChannelMarker[], SCBACollisionMarker[]
    fresh_scattering === nothing && return (; channels, collisions, status = :unavailable)
    mixed_scattering === nothing && return (; channels, collisions, status = :unavailable)
    for channel in problem.kernels.enabled
        mixed, fresh = mixed_scattering[channel], fresh_scattering[channel]
        _channel_component_markers!(channels, channel, mixed, fresh)
        push!(collisions, _channel_collision_marker(problem, green, channel, :mixed, mixed))
        push!(collisions, _channel_collision_marker(problem, green, channel, :fresh, fresh))
    end
    for (channel, mixed, fresh) in (
        (:embedding, mixed_embedding, fresh_embedding),
        (:embedding_plus, mixed_plus, fresh_plus),
        (:embedding_minus, mixed_minus, fresh_minus),
        (:total, total, fresh_total),
    )
        (mixed === nothing || fresh === nothing) && continue
        _channel_component_markers!(channels, channel, mixed, fresh)
    end
    return (; channels, collisions, status = :available)
end

"""Weight in edge bands whose ±shift requests leave the represented energy window.

This measures boundary exposure of the current state, not an extrapolated lost
collision flux. Signed weights are retained; PSD diagnostics assess validity.
"""
function _shift_boundary_markers(problem, green, shift)
    energy, weightsE, weightsK = problem.grids.ε, problem.grids.wᴱ, problem.grids.wᵏ
    occupied, spectral, edge_occupied, edge_spectral = 0.0, 0.0, 0.0, 0.0
    NE, Nk, Nb, _ = size(green.Gᴿ)
    @inbounds for e = 1:NE, m = 1:Nk
        weight = weightsE[e]*weightsK[m]
        n = sum(real(-im*green.Gˡ[e, m, a, a]) for a = 1:Nb)
        a = sum(real(green.A[e, m, j, j]) for j = 1:Nb)
        occupied += weight*n
        spectral += weight*a
        if energy[e]-abs(shift) < first(energy) || energy[e]+abs(shift) > last(energy)
            edge_occupied += weight*n
            edge_spectral += weight*a
        end
    end
    return (
        occupied > 0 ? edge_occupied/occupied : NaN,
        spectral > 0 ? edge_spectral/spectral : NaN,
    )
end

"""Reconstruct raw Keldysh blocks with O(Nb²) workspace; no full raw-array copy."""
function _scba_physical_markers(
    problem,
    green,
    total,
    current,
    iteration,
    policy::SCBAPhysicsMarkerPolicy,
    normalization::Symbol = :paired_convex;
    seed_mu_scaled::Float64 = NaN,
    seed_number_ratio::Float64 = NaN,
    map_context...,
)
    NE, Nk, Nb, _ = size(green.Gᴿ)
    grids = problem.grids
    target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    spin_factor = problem.physical.g_s/(2π)
    spectral_weights = zeros(NE)
    @inbounds for m = 1:Nk, e = 1:NE
        spectral_weights[e] +=
            spin_factor *
            grids.wᴱ[e] *
            grids.wᵏ[m] *
            sum(real(green.A[e, m, a, a]) for a = 1:Nb)
    end
    equilibrium_status = _marker_equilibrium_status(problem)
    equilibrium_applicable = equilibrium_status === :available
    mu = NaN
    kBT = _electronvolts(CODATA.kᴮₑᵥ*problem.physical.Tᴸ)/problem.scales.E₀_eV
    if equilibrium_applicable
        fit = _marker_spectral_chemical_potential(grids.ε, spectral_weights, target, kBT)
        mu, equilibrium_status = fit.mu, fit.status
    end
    GR, sigma, temporary, raw_n, raw_p = (zeros(ComplexF64, Nb, Nb) for _ = 1:5)
    raw_number = holes = norm_n = norm_p = correction_n = correction_p = 0.0
    fdt_raw_squared = fdt_normalized_squared = spectral_squared = 0.0
    seed_raw_squared = seed_normalized_squared = 0.0
    @inbounds for m = 1:Nk, e = 1:NE
        weight = spin_factor*grids.wᴱ[e]*grids.wᵏ[m]
        for b = 1:Nb, a = 1:Nb
            GR[a, b] = green.Gᴿ[e, m, a, b]
            sigma[a, b] = total.Σˡ[e, m, a, b]
        end
        mul!(temporary, GR, sigma)
        mul!(raw_n, temporary, adjoint(GR))
        for b = 1:Nb, a = 1:Nb
            sigma[a, b] = total.Σᵍ[e, m, a, b]
        end
        mul!(temporary, GR, sigma)
        mul!(raw_p, temporary, adjoint(GR))
        f = isfinite(mu) ? _fermi(grids.ε[e], mu, kBT) : NaN
        seed_f =
            equilibrium_applicable && isfinite(seed_mu_scaled) ?
            _fermi(grids.ε[e], seed_mu_scaled, kBT) : NaN
        for b = 1:Nb, a = 1:Nb
            gn, gp = -im*raw_n[a, b], im*raw_p[a, b]
            gn_normalized, gp_normalized = -im*green.Gˡ[e, m, a, b], im*green.Gᵍ[e, m, a, b]
            spectral = green.A[e, m, a, b]
            if a == b
                raw_number += weight*real(gn)
                holes += weight*real(gp)
            end
            norm_n += weight*abs2(gn)
            norm_p += weight*abs2(gp)
            correction_n += weight*abs2(gn_normalized-gn)
            correction_p += weight*abs2(gp_normalized-gp)
            spectral_squared += weight*abs2(spectral)
            if isfinite(seed_f)
                seed_raw_squared +=
                    weight*(abs2(gn-seed_f*spectral)+abs2(gp-(1-seed_f)*spectral))
                seed_normalized_squared +=
                    weight*(
                        abs2(gn_normalized-seed_f*spectral)+abs2(
                            gp_normalized-(1-seed_f)*spectral,
                        )
                    )
            end
            if isfinite(mu)
                fdt_raw_squared += weight*(abs2(gn-f*spectral)+abs2(gp-(1-f)*spectral))
                fdt_normalized_squared +=
                    weight*(
                        abs2(gn_normalized-f*spectral) + abs2(gp_normalized-(1-f)*spectral)
                    )
            end
        end
    end
    a, c = NaN, NaN
    if normalization === :paired_convex &&
       isfinite(raw_number) &&
       raw_number > 0 &&
       isfinite(holes) &&
       holes >= 0 &&
       target <= raw_number+holes
        coefficients = _paired_number_coefficients(raw_number, holes, target)
        a, c = coefficients.occupied_fraction, coefficients.empty_fraction
    elseif normalization === :scalar_lesser
        a, c = target/raw_number, 0.0
    end
    linewidth = _sampled_linewidth_markers(green.A, total, grids, policy)
    fresh = _fresh_map_markers(problem, green, total; map_context...)
    dE = (last(grids.ε)-first(grids.ε))/(NE-1)
    lo_shift = _electronvolts(problem.physical.ħωᴸᴼ)/problem.scales.E₀_eV
    field_shift = _scaled_physics(problem.physical, problem.scales).Eᵖ
    lo_edges =
        problem.scattering.LO ? _shift_boundary_markers(problem, green, lo_shift) :
        (NaN, NaN)
    field_edges = _shift_boundary_markers(problem, green, field_shift)
    return SCBAPhysicalMarkers(
        iteration,
        holes,
        raw_number+holes,
        a,
        c,
        _marker_relative_norm(correction_n, norm_n),
        _marker_relative_norm(correction_p, norm_p),
        equilibrium_applicable,
        equilibrium_status,
        mu*problem.scales.E₀_eV,
        isfinite(mu) ? _marker_relative_norm(fdt_raw_squared, spectral_squared) : NaN,
        isfinite(mu) ? _marker_relative_norm(fdt_normalized_squared, spectral_squared) :
        NaN,
        equilibrium_applicable ? abs(current) : NaN,
        (last(grids.ε)-first(grids.ε))/(NE-1)*problem.scales.E₀_eV,
        linewidth.q10,
        linewidth.q50,
        linewidth.q90,
        linewidth.underresolved,
        linewidth.blocks,
        linewidth.status,
        _MARKER_LINEWIDTH_SAMPLING_METHOD,
        policy.cadence,
        policy.max_spectral_blocks,
        policy.relative_mode_weight_floor,
        fresh.status,
        problem.scattering.LO ? lo_shift/dE : NaN,
        field_shift/dE,
        lo_edges...,
        field_edges...,
        seed_mu_scaled*problem.scales.E₀_eV,
        seed_number_ratio,
        equilibrium_applicable && isfinite(seed_mu_scaled) ?
        _marker_relative_norm(seed_raw_squared, spectral_squared) : NaN,
        equilibrium_applicable && isfinite(seed_mu_scaled) ?
        _marker_relative_norm(seed_normalized_squared, spectral_squared) : NaN,
        fresh.channels,
        fresh.collisions,
    )
end

function _with_scba_physical_markers(row::SCBAIteration, markers)
    return SCBAIteration(
        (
            getfield(row, name) for
            name in fieldnames(SCBAIteration) if name !== :physical_markers
        )...,
        markers,
    )
end

function _ensure_scba_physical_markers!(
    history,
    problem,
    green,
    total,
    options;
    terminal::Bool = false,
    map_context...,
)
    isempty(history) && return 0.0
    row = last(history)
    row.physical_markers !== nothing && return 0.0
    policy = options.physics_markers
    (terminal || row.ν == 1 || row.ν % policy.cadence == 0) || return 0.0
    started = time_ns()
    markers = _scba_physical_markers(
        problem,
        green,
        total,
        row.J,
        row.ν,
        policy,
        options.algorithms.occupation_normalization;
        map_context...,
    )
    history[end] = _with_scba_physical_markers(row, markers)
    return (time_ns()-started)*1e-9
end

"""First and sustained crossings on one ordered SCBA trajectory, never acceptance.

A missing iteration breaks a consecutive streak. Combining different outer
iterations or attempts by resetting iteration numbers is rejected. The owning
result supplies scientific/attempt/outer identity and final acceptance.
"""
function _scba_threshold_crossings(history; required_consecutive::Int = 3)
    required_consecutive > 0 || throw(ArgumentError("crossing streak must be positive"))
    all(i -> history[i].ν > history[i-1].ν, 2:length(history)) ||
        throw(ArgumentError("threshold crossings require one ordered SCBA trajectory"))
    summaries = NamedTuple[]
    value(row, metric) =
        metric === :joint ? max(row.r_Σ, row.r_K, row.r_λ) : getfield(row, metric)
    measurement(index, field) = iszero(index) ? NaN : getfield(history[index], field)
    iteration(index) = iszero(index) ? 0 : history[index].ν
    for metric in (:r_Σ, :r_K, :r_λ, :joint), threshold in (1e-3, 1e-4, 1e-5, 1e-6, 1e-8)
        first_index = sustained_start = sustained_end = streak = 0
        for index in eachindex(history)
            residual = value(history[index], metric)
            if isfinite(residual) && residual <= threshold
                iszero(first_index) && (first_index = index)
                contiguous = index > 1 && history[index].ν == history[index-1].ν+1
                streak = contiguous ? streak+1 : 1
                if streak >= required_consecutive && iszero(sustained_end)
                    sustained_start, sustained_end = index-required_consecutive+1, index
                end
            else
                streak = 0
            end
        end
        push!(
            summaries,
            (;
                metric,
                threshold,
                required_consecutive,
                first_iteration = iteration(first_index),
                sustained_start_iteration = iteration(sustained_start),
                sustained_end_iteration = iteration(sustained_end),
                first_current_A_m2 = measurement(first_index, :J),
                sustained_current_A_m2 = measurement(sustained_end, :J),
                first_raw_charge = measurement(first_index, :raw_charge),
                sustained_raw_charge = measurement(sustained_end, :raw_charge),
                first_population_change = measurement(first_index, :r_population),
                sustained_population_change = measurement(sustained_end, :r_population),
                first_PSD_residual = measurement(first_index, :r_PSD),
                sustained_PSD_residual = measurement(sustained_end, :r_PSD),
                first_causality_residual = measurement(first_index, :r_caus),
                sustained_causality_residual = measurement(sustained_end, :r_caus),
            ),
        )
    end
    return summaries
end
