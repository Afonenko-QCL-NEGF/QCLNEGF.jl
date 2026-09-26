function _electron_density_bar(problem::NEGFProblem, Gˡ::AbstractArray{<:Number,4})
    grids, χ = problem.grids, problem.basis.χ
    Nz, Nb = size(χ)
    size(Gˡ) == (length(grids.wᴱ), length(grids.wᵏ), Nb, Nb) ||
        throw(DimensionMismatch("lesser Green function and density basis disagree"))
    # Integrate every coherence once before projecting on the position grid.
    # No block copies or per-position E,k passes; no positivity repair.
    integrated = zeros(ComplexF64, Nb, Nb)
    for b = 1:Nb, a = 1:Nb
        value = zero(ComplexF64)
        for m in axes(Gˡ, 2), e in axes(Gˡ, 1)
            value += grids.wᴱ[e] * grids.wᵏ[m] * Gˡ[e, m, a, b]
        end
        integrated[a, b] = value
    end
    density = zeros(Float64, Nz)
    prefactor = problem.physical.g_s / (2π)
    for g = 1:Nz
        value = zero(ComplexF64)
        for b = 1:Nb, a = 1:Nb
            value += χ[g, a] * integrated[a, b] * conj(χ[g, b])
        end
        density[g] = prefactor * real(-im * value)
    end
    return density
end

"""
    electron_density(problem, Gˡ)

Reconstruct the physical volume density including all basis coherences.
Returns a Unitful vector in `m^-3`; internal Poisson code uses the private
dimensionless counterpart `L₀³n`.

Implements [EQ-DENSITY-001](@ref eq-realspace-density).
See [Green functions and density](@ref theory-greens) and
[observables and balances](@ref theory-observables).
"""
function electron_density(problem::NEGFProblem, Gˡ::AbstractArray{<:Number,4})
    return (_electron_density_bar(problem, Gˡ) ./ problem.scales.L₀_m^3) .* u"m^-3"
end

"""
    spatial_energy_density(problem, Gˡ)

Return `(energy_eV, z_m, n_per_eV_m3, axis_order)` with the density matrix
`n_per_eV_m3[energy, position]` in physical `eV⁻¹ m⁻³` and
`axis_order = (:energy, :position)`. Coordinates and densities are plain
floating-point arrays; their field names specify their units.

For the dimensionless Green array `(energy, momentum, basis, basis)`, the
projection is
`g_s/(2π E₀[eV] L₀[m]^3) * Σ_m wᵏ[m] Re[-i χ(z)ᵀ Gˡ(E,m) χ(z)*]`.
It retains every off-diagonal basis coherence and the radial momentum and
spin weights. Integrating the energy axis with `problem.grids.wᴱ .* E₀_eV`
recovers [`electron_density`](@ref); integrating the position axis with
`problem.grids.wˣ .* L₀_m` gives the energy-resolved sheet density.

The projection takes the real part of the density as in `electron_density`.
It neither clips negative values nor repairs an unphysical Green function.
`spectral_maps` retains its historical `(position, energy)` orientation and
uses the transpose of this same physical density.
"""
function spatial_energy_density(problem::NEGFProblem, Gˡ::AbstractArray{<:Number,4})
    grids = problem.grids
    χ = problem.basis.χ
    Nz, Nb = size(χ)
    NE, Nk = length(grids.ε), length(grids.wᵏ)
    size(Gˡ) == (NE, Nk, Nb, Nb) || throw(
        DimensionMismatch(
            "Gˡ must have shape (energy, momentum, basis, basis) = $((NE, Nk, Nb, Nb))",
        ),
    )
    length(grids.x) == Nz || throw(DimensionMismatch("basis and position grid disagree"))
    all(isfinite, Gˡ) || throw(ArgumentError("Gˡ must contain only finite values"))

    n = zeros(Float64, NE, Nz)
    integrated = zeros(ComplexF64, Nb, Nb)
    prefactor = problem.physical.g_s / (2π * problem.scales.E₀_eV * problem.scales.L₀_m^3)
    for e = 1:NE
        fill!(integrated, 0)
        for b = 1:Nb, a = 1:Nb, m = 1:Nk
            integrated[a, b] += grids.wᵏ[m] * Gˡ[e, m, a, b]
        end
        for z = 1:Nz
            value = 0.0 + 0.0im
            for b = 1:Nb, a = 1:Nb
                value += χ[z, a] * integrated[a, b] * conj(χ[z, b])
            end
            n[e, z] = prefactor * real(-im * value)
        end
    end
    return (
        energy_eV = grids.ε .* problem.scales.E₀_eV .+
                    _electronvolts(problem.physical.E_ref),
        z_m = grids.x .* problem.scales.L₀_m,
        n_per_eV_m3 = n,
        axis_order = (:energy, :position),
    )
end

function _sheet_density_matrix_bar(problem::NEGFProblem, Gˡ::AbstractArray{<:Number,4})
    Nb = size(Gˡ, 3)
    P = zeros(ComplexF64, Nb, Nb)
    for e in axes(Gˡ, 1), m in axes(Gˡ, 2)
        P .+=
            (-im * problem.physical.g_s / (2π)) *
            problem.grids.wᴱ[e] *
            problem.grids.wᵏ[m] .* _matrix_block(Gˡ, e, m)
    end
    return P
end

"""
Physical sheet-density matrix `P_ab` with Unitful `m^-2` entries.
Implements [EQ-POP-001](@ref eq-sheet-density).

See [Observables and balances](@ref theory-observables).
"""
function sheet_density_matrix(problem::NEGFProblem, Gˡ::AbstractArray{<:Number,4})
    return (_sheet_density_matrix_bar(problem, Gˡ) ./ problem.scales.L₀_m^2) .* u"m^-2"
end

"""
Return `(Nₛ,p)` for localized-state sheet populations and fractions.

See [Observables and balances](@ref theory-observables) and
[the sheet-density matrix](@ref eq-sheet-density).
"""
function state_populations(problem::NEGFProblem, Gˡ::AbstractArray{<:Number,4})
    Pbar = _sheet_density_matrix_bar(problem, Gˡ)
    Nbar = real.(diag(Pbar))
    target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    return (Nₛ = (Nbar ./ problem.scales.L₀_m^2) .* u"m^-2", p = Nbar ./ target)
end

"""
    effective_levels(problem, Uᴴ)

Diagonalize the full Hermitian projected Hamiltonian at `k=0` strictly as
post-processing.  Returns Unitful level energies and dipoles plus the
dimensionless one-band oscillator-strength matrix computed with
`problem.scales.m_ref`.  These diagnostics do not replace spectral peaks.

Implements [EQ-LEVELS-001](@ref eq-effective-levels).
See [Observables and balances](@ref theory-observables); these effective
levels are diagnostic eigenvalues, not interacting spectral peaks.
"""
function effective_levels(problem::NEGFProblem, Uᴴ::AbstractVector{<:Real})
    h₀ = Matrix(view(project_hamiltonians(problem, Uᴴ), 1, :, :))
    rH = norm(h₀ - h₀') / (norm(h₀) + _EFFECTIVE_HAMILTONIAN_NORM_FLOOR)
    rH ≤ _EFFECTIVE_HAMILTONIAN_HERMITICITY_LIMIT ||
        throw(DomainError(rH, "projected h(k=0) is not Hermitian"))
    decomposition = eigen(Hermitian(h₀))
    V = decomposition.vectors
    z_eff = V' * problem.basis.Z * V
    Nb = length(decomposition.values)
    f = zeros(Float64, Nb, Nb)
    mref = Float64(ustrip(u"kg", problem.scales.m_ref))
    ħ = Float64(ustrip(u"J*s", CODATA.ħ))
    for i = 1:Nb, j = 1:Nb
        ΔE = (decomposition.values[j] - decomposition.values[i]) * problem.scales.E₀_J
        zij = z_eff[i, j] * problem.scales.L₀_m
        f[i, j] = 2mref * ΔE * abs2(zij) / ħ^2
    end
    Eref = _electronvolts(problem.physical.E_ref)
    return (
        E = (decomposition.values .* problem.scales.E₀_eV .+ Eref) .* u"eV",
        z = z_eff .* problem.scales.L₀_m .* u"m",
        oscillator_strength = f,
        eigenvectors = V,
    )
end

function _current_trace(Σᵍ, Σˡ, Gˡ, Gᵍ, e, m)
    return tr(
        _matrix_block(Σᵍ, e, m) * _matrix_block(Gˡ, e, m) -
        _matrix_block(Σˡ, e, m) * _matrix_block(Gᵍ, e, m),
    )
end

function _energy_accumulation_order(problem::NEGFProblem; reverse_order::Bool = false)
    order = sortperm(
        collect(eachindex(problem.grids.ε));
        by = e -> (abs(problem.grids.ε[e]), problem.grids.ε[e]),
    )
    reverse_order && reverse!(order)
    return order
end

function _boundary_flux_complex(
    problem::NEGFProblem,
    family::SelfEnergyFamily,
    green::GreenState;
    reverse_order::Bool = false,
)
    value = 0.0 + 0.0im
    for e in _energy_accumulation_order(problem; reverse_order), m in axes(green.Gᴿ, 2)
        value +=
            problem.grids.wᴱ[e] *
            problem.grids.wᵏ[m] *
            _current_trace(family.Σᵍ, family.Σˡ, green.Gˡ, green.Gᵍ, e, m)
    end
    return problem.physical.g_s * value / (2π)
end

function _boundary_flux_bar(
    problem::NEGFProblem,
    family::SelfEnergyFamily,
    green::GreenState,
)
    return real(_boundary_flux_complex(problem, family, green))
end

"""
    boundary_flux(problem, family, green)

Outgoing electron number flux represented by one embedding family.  Returns
a Unitful value in `m^-2*s^-1`; positive `embedding_plus` means electron flow
toward `+z`.

Implements [EQ-CURRENT-001](@ref eq-boundary-flux).
See [Observables and balances](@ref theory-observables) for the current sign
and the relation between number flux and electrical current density.
"""
function boundary_flux(problem::NEGFProblem, family::SelfEnergyFamily, green::GreenState)
    flux₀ =
        problem.scales.E₀_eV /
        (Float64(ustrip(u"eV*s", CODATA.ħₑᵥ)) * problem.scales.L₀_m^2)
    return _boundary_flux_bar(problem, family, green) * flux₀ * u"m^-2*s^-1"
end

"""
Energy-resolved engineering current `j₊(E_e)` in `A*m^-2*eV^-1`.
Implements [EQ-CURRENT-002](@ref eq-resolved-current).

See [Observables and balances](@ref theory-observables).
"""
function energy_resolved_current(
    problem::NEGFProblem,
    family::SelfEnergyFamily,
    green::GreenState,
)
    NE, Nk = size(green.Gᴿ, 1), size(green.Gᴿ, 2)
    jbar = zeros(Float64, NE)
    for e = 1:NE, m = 1:Nk
        jbar[e] +=
            problem.grids.wᵏ[m] *
            real(_current_trace(family.Σᵍ, family.Σˡ, green.Gˡ, green.Gᵍ, e, m))
    end
    jbar .*= problem.physical.g_s / (2π)
    factor = Float64(ustrip(u"A/m^2", problem.scales.J₀)) / problem.scales.E₀_eV
    return jbar .* factor .* u"A/m^2/eV"
end

"""
    bond_current(problem, Gˡ)

Reconstruct the selected coordinate lesser blocks and return internal
electron-flow bond currents.  The periodic seam is intentionally absent and
must be compared with [`boundary_flux`](@ref).
Implements [EQ-CURRENT-002](@ref eq-resolved-current).
See [Observables and balances](@ref theory-observables) for its relation to
the embedding current and the omitted periodic seam.
"""
function bond_current(problem::NEGFProblem, Gˡ::AbstractArray{<:Number,4})
    Nz = problem.numerical.N_z
    Nb = problem.numerical.N_b
    Φ = problem.basis.Φ
    Hgrid = build_bdd_hamiltonian(problem.profiles, problem.grids, problem.scales)
    Jbar = zeros(Float64, Nz - 1)
    for g = 1:(Nz-1), e in axes(Gˡ, 1), m in axes(Gˡ, 2)
        Ggrid = 0.0 + 0.0im
        for a = 1:Nb, b = 1:Nb
            Ggrid += Φ[g+1, a] * Gˡ[e, m, a, b] * conj(Φ[g, b])
        end
        Jbar[g] += problem.grids.wᴱ[e] * problem.grids.wᵏ[m] * 2real(Hgrid[g, g+1] * Ggrid)
    end
    Jbar .*= problem.physical.g_s / (2π)
    return Jbar .* Float64(ustrip(u"A/m^2", problem.scales.J₀)) .* u"A/m^2"
end

function _collision_sum(
    problem::NEGFProblem,
    family::SelfEnergyFamily,
    green::GreenState;
    reverse_order::Bool = false,
)
    signed = 0.0 + 0.0im
    absolute = 0.0
    for e in _energy_accumulation_order(problem; reverse_order), m in axes(green.Gᴿ, 2)
        incoming = tr(_matrix_block(family.Σˡ, e, m) * _matrix_block(green.Gᵍ, e, m))
        outgoing = tr(_matrix_block(family.Σᵍ, e, m) * _matrix_block(green.Gˡ, e, m))
        weight = problem.grids.wᴱ[e] * problem.grids.wᵏ[m]
        signed += weight * (incoming - outgoing)
        absolute += weight * (abs(incoming) + abs(outgoing))
    end
    return signed, absolute
end

"""
Particle-conservation rate and normalized residual for one mechanism.
The returned diagnostics also expose the discarded imaginary part and the
forward/reverse fixed-order cancellation estimate.
Implements [EQ-COLLISION-001](@ref eq-collision-balance).

See [Observables and balances](@ref theory-observables).
"""
function collision_balance(
    problem::NEGFProblem,
    family::SelfEnergyFamily,
    green::GreenState,
)
    signed, absolute = _collision_sum(problem, family, green)
    signed_reverse, _ = _collision_sum(problem, family, green; reverse_order = true)
    prefbar = problem.physical.g_s / (2π)
    floor = 1e-14
    residual = prefbar * abs(real(signed)) / (prefbar * absolute + floor)
    roundoff =
        abs(signed) > 1e-12 ? abs(signed - signed_reverse) / abs(signed) :
        100 * abs(signed - signed_reverse)
    flux₀ =
        problem.scales.E₀_eV /
        (Float64(ustrip(u"eV*s", CODATA.ħₑᵥ)) * problem.scales.L₀_m^2)
    return (
        rate = prefbar * real(signed) * flux₀ * u"m^-2*s^-1",
        residual = residual,
        imaginary = imag(signed),
        roundoff = roundoff,
    )
end

function _collision_power_bar(
    problem::NEGFProblem,
    family::SelfEnergyFamily,
    green::GreenState;
    reverse_order::Bool = false,
)
    value = 0.0 + 0.0im
    for e in _energy_accumulation_order(problem; reverse_order), m in axes(green.Gᴿ, 2)
        incoming = tr(_matrix_block(family.Σˡ, e, m) * _matrix_block(green.Gᵍ, e, m))
        outgoing = tr(_matrix_block(family.Σᵍ, e, m) * _matrix_block(green.Gˡ, e, m))
        value +=
            problem.grids.wᴱ[e] *
            problem.grids.wᵏ[m] *
            problem.grids.ε[e] *
            (incoming - outgoing)
    end
    return problem.physical.g_s * value / (2π)
end

_power_scale(problem::NEGFProblem) =
    problem.scales.E₀_J * problem.scales.E₀_eV /
    (Float64(ustrip(u"eV*s", CODATA.ħₑᵥ)) * problem.scales.L₀_m^2)

"""Energy transferred from one mechanism to electrons, in `W/m²`."""
function _collision_power(problem::NEGFProblem, family::SelfEnergyFamily, green::GreenState)
    value = _collision_power_bar(problem, family, green)
    return real(value) * _power_scale(problem) * u"W/m^2"
end

"""
Total electronic power-balance residual over internal scattering and bias
work. Implements [EQ-POWER-001](@ref eq-power-balance).

See [Observables and balances](@ref theory-observables).
"""
function power_balance(
    problem::NEGFProblem,
    mechanisms::Dict{Symbol,SelfEnergyFamily},
    embedding_plus::SelfEnergyFamily,
    green::GreenState,
)
    powers = Dict(
        name => _collision_power(problem, family, green) for (name, family) in mechanisms
    )
    Psc = sum(Float64(ustrip(u"W/m^2", P)) for P in values(powers); init = 0.0)
    J = CODATA.e * boundary_flux(problem, embedding_plus, green)
    Vp = uconvert(u"V", problem.physical.F_bias * period_length(problem.physical))
    Pfield = Float64(ustrip(u"W/m^2", uconvert(u"W/m^2", J * Vp)))
    P₀ = _power_scale(problem)
    residual =
        abs(Psc + Pfield) / (
            sum(abs(Float64(ustrip(u"W/m^2", P))) for P in values(powers); init = 0.0) +
            abs(Pfield) +
            1e-14 * P₀
        )
    imaginary = 0.0
    roundoff = 0.0
    for family in values(mechanisms)
        forward = _collision_power_bar(problem, family, green)
        reversed = _collision_power_bar(problem, family, green; reverse_order = true)
        imaginary = max(imaginary, abs(imag(forward)))
        estimate =
            abs(forward) > 1e-12 ? abs(forward - reversed) / abs(forward) :
            100 * abs(forward - reversed)
        roundoff = max(roundoff, estimate)
    end
    return (
        mechanisms = powers,
        field = Pfield * u"W/m^2",
        residual = residual,
        imaginary = imaginary,
        roundoff = roundoff,
    )
end

"""
Return physical local spectral and occupied maps for visualization.

See [Green functions and density](@ref theory-greens) for their definitions
and [debug/final visualization](@ref theory-visualization) for their use.
"""
function spectral_maps(problem::NEGFProblem, green::GreenState)
    Nz = problem.numerical.N_z
    NE = problem.numerical.N_E
    χ = problem.basis.χ
    Amap = zeros(Float64, Nz, NE)
    density = spatial_energy_density(problem, green.Gˡ)
    for g = 1:Nz, e = 1:NE, m = 1:problem.numerical.N_k
        χg = view(χ, g, :)
        Amap[g, e] +=
            problem.physical.g_s *
            problem.grids.wᵏ[m] *
            real(transpose(χg) * _matrix_block(green.A, e, m) * conj.(χg))
    end
    density_unit = 1 / (problem.scales.E₀_eV * problem.scales.L₀_m^3)
    return (
        z = problem.grids.x .* problem.scales.L₀_m .* u"m",
        E = (
            problem.grids.ε .* problem.scales.E₀_eV .+
            _electronvolts(problem.physical.E_ref)
        ) .* u"eV",
        A = Amap .* density_unit .* u"eV^-1*m^-3",
        N = permutedims(density.n_per_eV_m3) .* u"eV^-1*m^-3",
    )
end
