# Representative-momentum SPPA correlation operator. The same physical
# approximation is available to reference and optimized SCBA. Physical limits
# are tested independently; agreement of these two routes is not an oracle.

function sppa_parameters(physical, options)
    options.mode === :sppa_single_q || throw(ArgumentError("SPPA model is not enabled"))
    physical.g_s == 2 || throw(
        ArgumentError("SPPA implementation assumes spin-degenerate unpolarized electrons"),
    )
    mass = options.effective_mass_ratio * CODATA.m₀
    screening = if options.screening_dimension == 2
        PhysicalModelExtensions.ThomasFermiScreening2D(
            options.carrier_density_si * u"m^-2",
            options.electron_temperature,
            mass,
            physical.g_s,
            physical.ε_s,
        )
    else
        PhysicalModelExtensions.DebyeScreening3D(
            options.carrier_density_si * u"m^-3",
            options.electron_temperature,
            physical.ε_s,
        )
    end
    return PhysicalModelExtensions.single_plasmon_pole(
        screening,
        options.transfer_wavenumber,
        mass,
        options.electron_temperature,
    )
end

"""
Construct the full four-index quasi-2D subband form factor at a declared q.
The two spatial quadratures use the existing basis normalization. Neither
diagonal-state reduction nor empirical relaxation times are introduced.
"""
function sppa_form_factor(grids::ModelGrids, basis::BasisData, transfer_dimensionless::Real)
    transfer_dimensionless > 0 && isfinite(transfer_dimensionless) ||
        throw(ArgumentError("representative transfer momentum must be positive"))
    nz, nb = size(basis.χ)
    length(grids.x) == nz && length(grids.wˣ) == nz ||
        throw(DimensionMismatch("basis/grid axes differ"))
    vertices = Matrix{ComplexF64}(undef, nz, nb^2)
    for c = 1:nb, a = 1:nb, g = 1:nz
        vertices[g, a+(c-1)*nb] = conj(basis.χ[g, a]) * basis.χ[g, c] * grids.wˣ[g]
    end
    metric = [exp(-transfer_dimensionless * abs(x-y)) for x in grids.x, y in grids.x]
    # V_ac V_bd* ordering creates a positive covariance for the map K_acdb G_cd.
    covariance = transpose(vertices) * metric * conj(vertices)
    form = zeros(ComplexF64, nb, nb, nb, nb)
    for b = 1:nb, d = 1:nb, c = 1:nb, a = 1:nb
        form[a, c, d, b] = covariance[a+(c-1)*nb, b+(d-1)*nb]
    end
    return form
end

"""
Add the actual normalized six-axis SPPA kernel to KernelSet. The bare Coulomb
factor remains quasi-2D for either screening closure, following Winge (5)–(9).
The declared representative q replaces angular/momentum-transfer dependence.
"""
function add_sppa_kernel(
    physical::PhysicalParameters,
    scales::ScaleSystem,
    grids::ModelGrids,
    basis::BasisData,
    kernels::KernelSet,
    options,
)
    options.mode === :none && return kernels
    :electron_electron in kernels.enabled &&
        throw(ArgumentError("e-e kernel is already installed"))
    pole = sppa_parameters(physical, options)
    q = _inverse_metres(options.transfer_wavenumber)
    form = sppa_form_factor(grids, basis, q * scales.L₀_m)
    bare = uconvert(
        u"J*m^2",
        CODATA.e^2 / (2 * CODATA.ε₀ * physical.ε_s * options.transfer_wavenumber),
    )
    factor = Float64(ustrip(u"J^2*m^2", bare * pole.correlation_weight)) / _K₀(scales)
    finite_form = factor .* form
    normalization = maximum(abs, finite_form)
    isfinite(normalization) && normalization > 0 ||
        throw(ArgumentError("SPPA kernel normalization must be finite and positive"))
    nk, nb = length(grids.κ), size(basis.χ, 2)
    tensor = Array{ComplexF64,6}(undef, nk, nk, nb, nb, nb, nb)
    normalized = finite_form ./ normalization
    for m = 1:nk, mp = 1:nk
        tensor[m, mp, :, :, :, :] .= normalized
    end
    definitions, factors = copy(kernels.K), copy(kernels.qᴷ)
    definitions[:electron_electron] = tensor
    factors[:electron_electron] = normalization
    return KernelSet(
        definitions,
        factors,
        kernels.Fᴸᴼ,
        vcat(kernels.enabled, :electron_electron),
    )
end

function _sppa_local_contraction(problem::NEGFProblem, field::Array{ComplexF64,4})
    ne, nk, nb, nb2 = size(field)
    nb == nb2 || throw(DimensionMismatch("Green blocks must be square"))
    weights = problem.grids.wᵏ
    length(weights) == nk || throw(DimensionMismatch("momentum weights differ"))
    averaged = zeros(ComplexF64, ne, nb, nb)
    for d = 1:nb, c = 1:nb, m = 1:nk, e = 1:ne
        averaged[e, c, d] += weights[m] * field[e, m, c, d]
    end
    result = zeros(ComplexF64, ne, nb, nb)
    kernel = view(problem.kernels.K[:electron_electron], 1, 1, :, :, :, :)
    scale = problem.kernels.qᴷ[:electron_electron]
    for b = 1:nb, a = 1:nb, d = 1:nb, c = 1:nb, e = 1:ne
        result[e, a, b] += scale * kernel[a, c, d, b] * averaged[e, c, d]
    end
    return result
end

"""Bounded positive shift pair without dense NE×NE matrices or cyclic wrapping."""
function _sppa_shift_pair(
    grids::ModelGrids,
    displacement::Real,
    field::Array{ComplexF64,3},
    discretization::Symbol,
)
    energies, weights = grids.ε, grids.wᴱ
    ne, nb, _ = size(field)
    length(energies) == ne || throw(DimensionMismatch("energy axes differ"))
    plus, minus = zeros(ComplexF64, size(field)), zeros(ComplexF64, size(field))
    if discretization === :finite_volume_piecewise_constant
        edges =
            vcat(energies[1], (energies[1:(end-1)] .+ energies[2:end]) ./ 2, energies[end])
        for i = 1:ne
            left, right = edges[i] + displacement, edges[i+1] + displacement
            j = max(1, searchsortedlast(edges, left))
            while j <= ne && edges[j] < right
                overlap = max(0.0, min(right, edges[j+1]) - max(left, edges[j]))
                if overlap > 0
                    @views plus[i, :, :] .+= (overlap / weights[i]) .* field[j, :, :]
                    @views minus[j, :, :] .+= (overlap / weights[j]) .* field[i, :, :]
                end
                j += 1
            end
        end
    elseif discretization === :nodal_linear
        for (offset, shifted) in ((displacement, plus), (-displacement, minus)), i = 1:ne
            target = energies[i] + offset
            energies[1] <= target <= energies[end] || continue
            j = clamp(searchsortedlast(energies, target), 1, ne-1)
            fraction = (target - energies[j]) / (energies[j+1] - energies[j])
            @views shifted[i, :, :] .=
                (1-fraction) .* field[j, :, :] .+ fraction .* field[j+1, :, :]
        end
    else
        throw(ArgumentError("unsupported SPPA energy-shift discretization"))
    end
    return plus, minus
end

"""
Return dimensionless lesser, greater and exchange components for actual SCBA.
The dynamic Hilbert transform is left to the caller's selected numerical
backend. Exchange is the bare instantaneous Fock term and is added only to ΣR.
"""
function sppa_components(problem::NEGFProblem, green::GreenState)
    options = problem.models.electron_electron
    pole = sppa_parameters(problem.physical, options)
    :electron_electron in problem.kernels.enabled ||
        throw(ArgumentError("SPPA kernel has not been assembled"))
    less = _sppa_local_contraction(problem, green.Gˡ)
    greater = _sppa_local_contraction(problem, green.Gᵍ)
    delta = _electronvolts(pole.pole_energy) / problem.scales.E₀_eV
    lp, lm =
        _sppa_shift_pair(problem.grids, delta, less, problem.energy_shift_discretization)
    gp, gm =
        _sppa_shift_pair(problem.grids, delta, greater, problem.energy_shift_discretization)
    n = pole.occupation
    local_lesser = (n+1) .* lp .+ n .* lm
    local_greater = (n+1) .* gm .+ n .* gp
    sigma_lesser = zeros(ComplexF64, size(green.Gˡ))
    sigma_greater = similar(sigma_lesser)
    for m in axes(sigma_lesser, 2)
        sigma_lesser[:, m, :, :] .= local_lesser
        sigma_greater[:, m, :, :] .= local_greater
    end
    nb = size(green.Gˡ, 3)
    exchange = zeros(ComplexF64, nb, nb)
    if options.include_exchange
        dynamic_energy = _electronvolts(pole.correlation_weight) / problem.scales.E₀_eV
        for b = 1:nb, a = 1:nb, e in eachindex(problem.grids.ε)
            exchange[a, b] +=
                im / (2π * dynamic_energy) * problem.grids.wᴱ[e] * less[e, a, b]
        end
        mismatch = norm(exchange - exchange') / max(norm(exchange), floatmin(Float64))
        mismatch <= 1e-9 || throw(
            DomainError(
                mismatch,
                "SPPA Fock shift is not Hermitian; input density is inconsistent",
            ),
        )
        exchange = (exchange + exchange') / 2
    end
    return (
        lesser = sigma_lesser,
        greater = sigma_greater,
        exchange = exchange,
        pole = pole,
    )
end

"""Reference Hilbert completion of the SPPA family, including the Fock term."""
function sppa_self_energy(
    problem::NEGFProblem,
    green::GreenState;
    return_roundoff::Bool = false,
)
    components = sppa_components(problem, green)
    if return_roundoff
        retarded, _, roundoff = retarded_self_energy(
            components.lesser,
            components.greater,
            problem.grids;
            return_roundoff = true,
        )
    else
        retarded, _ =
            retarded_self_energy(components.lesser, components.greater, problem.grids)
        roundoff = 0.0
    end
    for m in axes(retarded, 2), e in axes(retarded, 1)
        @views retarded[e, m, :, :] .+= components.exchange
    end
    family = SelfEnergyFamily(retarded, components.lesser, components.greater)
    return return_roundoff ? (family, roundoff) : family
end

"""
Measure charge and electronic-energy exchange of an e-e family. An equilibrium
plasmon closure may exchange energy with that effective bath; a small SCBA
residual does not certify this as conserving electronic GW transport.
"""
function electron_electron_collision_diagnostics(
    problem::NEGFProblem,
    green::GreenState,
    family::SelfEnergyFamily,
)
    number, energy, number_scale, energy_scale = 0.0, 0.0, 0.0, 0.0
    for e in eachindex(problem.grids.ε), m in eachindex(problem.grids.κ)
        incoming = real(tr(_matrix_block(family.Σˡ, e, m) * _matrix_block(green.Gᵍ, e, m)))
        outgoing = real(tr(_matrix_block(family.Σᵍ, e, m) * _matrix_block(green.Gˡ, e, m)))
        weight = problem.grids.wᴱ[e] * problem.grids.wᵏ[m]
        collision = incoming - outgoing
        activity = abs(incoming) + abs(outgoing)
        number += weight * collision
        energy += weight * problem.grids.ε[e] * collision
        number_scale += weight * activity
        energy_scale += weight * abs(problem.grids.ε[e]) * activity
    end
    number_factor =
        problem.physical.g_s * problem.scales.E₀_J /
        (2π * Float64(ustrip(u"J*s", CODATA.ħ)) * problem.scales.L₀_m^2)
    return (
        charge_residual = abs(number) / max(number_scale, floatmin(Float64)),
        energy_residual = abs(energy) / max(energy_scale, floatmin(Float64)),
        electron_number_rate_m2_s = number_factor * number,
        electron_energy_rate_W_m2 = number_factor * problem.scales.E₀_J * energy,
        closure = :equilibrium_plasmon_at_declared_electron_temperature,
        conserving_GW_certificate = false,
    )
end
