#!/usr/bin/env julia
# Reproduce the documentation data with the real domain/numerics layer only.
# julia --startup-file=no --project=. examples/10_physics_atlas.jl [output-directory]
using QCLNEGF
using LinearAlgebra
using Unitful
using Printf

const OUTPUT = abspath(
    isempty(ARGS) ?
    joinpath(@__DIR__, "..", "docs", "src", "assets", "physics", "data") : only(ARGS),
)
mkpath(OUTPUT)
BLAS.set_num_threads(1)

function csv(name, header, rows)
    open(joinpath(OUTPUT, name), "w") do io
        println(io, join(header, ','))
        for row in rows
            println(
                io,
                join(
                    (x isa AbstractFloat ? @sprintf("%.12g", x) : string(x) for x in row),
                    ',',
                ),
            )
        end
    end
end

function numerics(Nz = 96, P = 2, Nb = 5)
    NumericalParameters(
        N_z = Nz,
        N_b = Nb,
        P_basis = P,
        E_min = -0.2u"eV",
        E_max = 0.5u"eV",
        N_E = 31,
        M_E = 0.08u"eV",
        k_max = 0.4u"nm^-1",
        N_k = 3,
        N_φ = 8,
        qz_max = 5u"nm^-1",
        N_qz = 9,
        η_seed = 2u"meV",
    )
end

p, scales, n = reference_parameters(), ScaleSystem(), numerics()
grids = build_grids(p, n, scales)
profiles = build_profiles(p, n, scales, grids)
sp = QCLNEGF._scaled_physics(p, scales)
nm = scales.L₀_m * 1e9
meV = scales.E₀_eV * 1e3
dx_nm = grids.wˣ[1] * nm
period_nm = sp.Lp * nm
hbar_meV_ps = ustrip(u"meV*ps", QCLNEGF.CODATA.ħₑᵥ)
layer_rows = Tuple[]
let left = 0.0
    for (index, layer) in enumerate(p.layers)
        right = left + ustrip(u"nm", layer.d)
        push!(layer_rows, (index, left, right, layer.material, layer.doped))
        left = right
    end
end
csv("layers.csv", ["layer", "left_nm", "right_nm", "material", "doped"], layer_rows)
csv(
    "parameters.csv",
    ["name", "value", "unit"],
    [
        ("period", period_nm, "nm"),
        ("spatial_points_per_period", n.N_z, "1"),
        ("basis_states_per_period", n.N_b, "1"),
        ("periods", 2*n.P_basis+1, "1"),
        ("phonon_energy", ustrip(u"meV", p.ħωᴸᴼ), "meV"),
        ("reduced_planck_constant", hbar_meV_ps, "meV*ps"),
        ("boltzmann_constant", ustrip(u"meV/K", QCLNEGF.CODATA.kᴮₑᵥ), "meV/K"),
    ],
)
csv(
    "structure.csv",
    ["z_nm", "Ec_meV", "bias_meV", "mass_ratio", "donors_m3"],
    (
        (
            grids.x[g]*nm,
            profiles.Eᶜ[g]*meV,
            -sp.F*grids.x[g]*meV,
            profiles.mᶻᵣ[g],
            profiles.Nᴰ[g]/scales.L₀_m^3,
        ) for g in eachindex(grids.x)
    ),
)

references = Dict(
    method => QCLNEGF.build_multiperiod_reference_basis(
        p,
        n,
        scales,
        grids,
        profiles;
        localization = method,
    ) for method in (:pzp_tails, :bloch_wannier, :wannier_stark)
)
reference = references[:pzp_tails]
eig = eigen(Hermitian(reference.real_space_hamiltonian))
Nstates = length(reference.reference_energies)
C = ComplexF64.(eig.vectors[:, 1:Nstates])
QCLNEGF._phase_fix!(C)
energies = eig.values[1:Nstates] .* meV

basis_rows = Tuple[]
summary_rows = Tuple[]
operator_rows = Tuple[]
for method in (:energy, :pzp_tails, :bloch_wannier, :wannier_stark)
    ref = method === :energy ? reference : references[method]
    transform = method === :energy ? C : ref.transform
    H = transform' * ref.real_space_hamiltonian * transform
    Z = transform' * Diagonal(ref.coordinates) * transform
    v = im .* (H * Z - Z * H) .* (meV * nm / hbar_meV_ps) # nm/ps
    for j in axes(transform, 2)
        centre = real(Z[j, j]) * nm
        spread = sqrt(
            max(sum(abs2.(transform[:, j]) .* (ref.coordinates .* nm .- centre) .^ 2), 0.0),
        )
        ipr = sum(abs2.(transform[:, j]) .^ 2) / dx_nm
        tail = 1 - sum(abs2, transform[ref.central_rows, j])
        push!(
            summary_rows,
            (
                method,
                j,
                centre,
                spread,
                ipr,
                tail,
                real(H[j, j])*meV,
                norm(transform'*transform-I),
                norm(ref.real_space_hamiltonian*transform-transform*H)/norm(
                    ref.real_space_hamiltonian*transform,
                ),
            ),
        )
        # Keep the five displayed envelopes, while metrics/operators retain
        # every selected state. Full arrays need not bloat documentation.
        if j in ref.central_columns
            for g in eachindex(ref.coordinates)
                push!(
                    basis_rows,
                    (
                        method,
                        j,
                        ref.coordinates[g]*nm,
                        real(transform[g, j])/sqrt(dx_nm),
                        imag(transform[g, j])/sqrt(dx_nm),
                        abs2(transform[g, j])/dx_nm,
                    ),
                )
            end
        end
    end
    for a in axes(H, 1), b in axes(H, 2)
        push!(
            operator_rows,
            (
                method,
                a,
                b,
                real(H[a, b])*meV,
                imag(H[a, b])*meV,
                real(Z[a, b])*nm,
                imag(Z[a, b])*nm,
                abs(v[a, b]),
            ),
        )
    end
    @assert norm(transform'*transform-I) < 1e-10
    @assert norm(ref.real_space_hamiltonian*transform-transform*H)/norm(
        ref.real_space_hamiltonian*transform,
    ) < 1e-10
end
csv(
    "basis_envelopes.csv",
    [
        "method",
        "state",
        "z_nm",
        "real_psi_per_sqrt_nm",
        "imag_psi_per_sqrt_nm",
        "probability_per_nm",
    ],
    basis_rows,
)
csv(
    "basis_metrics.csv",
    [
        "method",
        "state",
        "centre_nm",
        "spread_nm",
        "ipr_per_nm",
        "outside_central_weight",
        "mean_energy_meV",
        "orthogonality",
        "subspace_residual",
    ],
    summary_rows,
)
csv(
    "basis_operators.csv",
    [
        "method",
        "a",
        "b",
        "H_real_meV",
        "H_imag_meV",
        "Z_real_nm",
        "Z_imag_nm",
        "speed_nm_ps",
    ],
    operator_rows,
)

# Same finite-dimensional resolvent in two bases. This is a noninteracting
# spectral illustration: eta and the Fermi distribution are prescribed.
eta = ustrip(u"meV", 2.0u"meV")
chemical_potential = ustrip(u"meV", 60.0u"meV")
kBT = ustrip(u"meV", QCLNEGF.CODATA.kᴮₑᵥ*p.Tᴸ)
# Explicit educational probe energy, not a solver default.
probe = ustrip(u"meV", 73.0u"meV")
H_meV = reference.hamiltonian .* meV
G_projected = inv((probe + im*eta)*I - H_meV)
G_coordinate = reference.transform * G_projected * reference.transform'
G_energy = C * Diagonal(1 ./ (probe .+ im*eta .- energies)) * C'
resolvent_error = norm(G_coordinate-G_energy)/norm(G_energy)
@assert resolvent_error < 1e-10
energy_grid = range(-20.0, 230.0; length = 181)
spectral_rows = Tuple[]
for E in energy_grid
    lorentz = eta ./ (π .* ((E .- energies) .^ 2 .+ eta^2))
    f = 1/(1+exp((E-chemical_potential)/kBT))
    density_of_states = abs2.(C) * lorentz / dx_nm
    for g in reference.central_rows
        push!(
            spectral_rows,
            (E, reference.coordinates[g]*nm, density_of_states[g], f*density_of_states[g]),
        )
    end
end
csv(
    "spectral_population.csv",
    ["energy_meV", "z_nm", "ldos_per_meV_nm", "occupied_per_meV_nm"],
    spectral_rows,
)

# A neutral trial electron profile constructed from actual localized states.
# One Poisson solve measures electrostatic response, not self-consistency.
central = reference.transform[reference.central_rows, reference.central_columns]
trial_shape = abs2.(central) * [0.4, 0.3, 0.15, 0.10, 0.05]
electron_density = sp.Nᴰ²ᴰ .* trial_shape ./ dot(grids.wˣ, trial_shape)
UH, gauge_multiplier, poisson_residual, neutrality =
    QCLNEGF.solve_periodic_poisson(profiles, grids, scales, electron_density)
@assert poisson_residual < 1e-10 && neutrality < 1e-12
@assert abs(sum(UH)) < 1e-10
csv(
    "hartree.csv",
    ["z_nm", "donors_m3", "electrons_m3", "hartree_meV", "total_band_meV"],
    (
        (
            grids.x[g]*nm,
            profiles.Nᴰ[g]/scales.L₀_m^3,
            electron_density[g]/scales.L₀_m^3,
            UH[g]*meV,
            (profiles.Eᶜ[g]-sp.F*grids.x[g]+UH[g])*meV,
        ) for g in eachindex(grids.x)
    ),
)

# Longitudinal form factors of full-window PzP states: exact projected exp(iqz).
a, b = first(reference.central_columns), first(reference.central_columns)+1
screening_q_nm = ustrip(u"nm^-1", p.q_s)
correlation_nm = ustrip(u"nm", p.Λᴵᶠᴿ)
csv(
    "form_factors.csv",
    [
        "q_per_nm",
        "diagonal_F2",
        "transition_F2",
        "normalized_screening",
        "normalized_gaussian_IFR",
    ],
    (
        (
            q,
            abs2(
                dot(
                    reference.transform[:, a],
                    exp.(im*q .* reference.coordinates .* nm) .* reference.transform[:, a],
                ),
            ),
            abs2(
                dot(
                    reference.transform[:, a],
                    exp.(im*q .* reference.coordinates .* nm) .* reference.transform[:, b],
                ),
            ),
            screening_q_nm^2/(q^2+screening_q_nm^2),
            exp(-q^2*correlation_nm^2/4),
        ) for q in range(0.0, 1.2; length = 241)
    ),
)

# Grid/window diagnostics compare invariant subspaces, not transport currents.
refinement_rows = Tuple[]
for (Nz, P) in ((48, 2), (96, 2), (192, 2), (96, 1), (96, 3))
    ni = numerics(Nz, P)
    gi = build_grids(p, ni, scales)
    pi = build_profiles(p, ni, scales, gi)
    ri = QCLNEGF.build_multiperiod_reference_basis(
        p,
        ni,
        scales,
        gi,
        pi;
        localization = :pzp_tails,
    )
    for (j, c) in enumerate(ri.central_columns)
        push!(
            refinement_rows,
            (
                Nz,
                2P+1,
                j,
                ri.centres[c]*nm,
                ri.spreads[c]*nm,
                ri.central_tail_weights[j],
                ri.subspace_residual,
                maximum(abs, eigvals(Hermitian(ri.hamiltonian))-ri.reference_energies)*meV,
            ),
        )
    end
end
csv(
    "localization_refinement.csv",
    [
        "Nz_per_period",
        "periods",
        "central_state",
        "centre_nm",
        "spread_nm",
        "tail_weight",
        "subspace_residual",
        "spectral_error_meV",
    ],
    refinement_rows,
)

# A two-level line shape uses a real cell Hamiltonian/dipole but prescribed
# occupations and linewidth. It is neither a cavity spectrum nor QCL gain.
cell_H = build_bdd_hamiltonian(profiles, grids, scales)
cell_eig = eigen(Hermitian(cell_H))
cell_Z = cell_eig.vectors' * Diagonal(grids.x .* nm) * cell_eig.vectors
transition_meV = (cell_eig.values[2]-cell_eig.values[1])*meV
dipole_nm = abs(cell_Z[1, 2])
csv(
    "optical_parameters.csv",
    [
        "transition_meV",
        "dipole_nm",
        "halfwidth_meV",
        "temperature_K",
        "spectral_eta_meV",
        "chemical_potential_meV",
    ],
    [(transition_meV, dipole_nm, 2.0, 200.0, eta, chemical_potential)],
)
csv(
    "dispersion.csv",
    ["k_per_nm", "state", "energy_meV"],
    (
        (k, j, energy*meV) for k in range(0.0, 0.4; length = 81) for
        (j, energy) in enumerate(
            eigvals(
                Hermitian(
                    cell_H +
                    Diagonal(QCLNEGF._C̄(scales)*(k*nm)^2 ./ profiles.m_parallelᵣ),
                ),
            )[1:5],
        )
    ),
)
csv(
    "checks.csv",
    ["check", "value", "tolerance"],
    [
        ("resolvent_basis_invariance", resolvent_error, 1e-10),
        ("poisson_residual", poisson_residual, 1e-10),
        ("charge_neutrality", neutrality, 1e-12),
        ("hartree_mean_dimensionless", abs(sum(UH)/length(UH)), 1e-10),
    ],
)
println("Physics atlas data: ", OUTPUT)
println(
    "Actual reference design profiles and finite-window operators; prescribed broadening/occupations; no converged transport run.",
)
