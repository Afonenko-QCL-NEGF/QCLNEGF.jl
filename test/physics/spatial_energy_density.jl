module SpatialEnergyDensityTests
include("../support/common.jl")
using QCLNEGF: spatial_energy_density

function replace_fields(value, changes::NamedTuple)
    return typeof(value)(
        (
            get(changes, field, getfield(value, field)) for
            field in fieldnames(typeof(value))
        )...,
    )
end

@testset "Physical spatial-energy density: units, coherence and integrated populations" begin
    scattering = ScatteringOptions(
        LO = false,
        acoustic = false,
        impurity = false,
        IFR = false,
        alloy = false,
    )
    original = build_problem(
        physical = reference_parameters(Tᴸ = 70u"K"),
        numerical = tutorial_numerics(),
        scattering = scattering,
        validate_static = false,
    )
    ne, nk, nb = original.numerical.N_E, original.numerical.N_k, original.numerical.N_b
    nz = length(original.grids.x)

    # Arbitrary complex rephasing is a physical representation change. It
    # makes a mistaken conjugation of a basis function observable here.
    phases = Diagonal(cis.(range(0.2, 1.6; length = nb)))
    basis = original.basis
    phased_basis = replace_fields(
        basis,
        (
            Φ = basis.Φ*phases,
            χ = basis.χ*phases,
            H₀ = phases'*basis.H₀*phases,
            Z = phases'*basis.Z*phases,
            M⁻¹ = phases'*basis.M⁻¹*phases,
            T₊ = phases'*basis.T₊*phases,
            T₋ = phases'*basis.T₋*phases,
        ),
    )
    problem = replace_fields(original, (basis = phased_basis,))
    v = normalize(ComplexF64[1 + 0.3im*a for a = 1:nb])
    phased_v = phases' * v
    energy_amplitude = [0.4 + exp(-0.01*(e-ne/2)^2) for e = 1:ne]
    momentum_amplitude = [0.7 + 0.2m for m = 1:nk]
    lesser = zeros(ComplexF64, ne, nk, nb, nb)
    original_lesser = similar(lesser)
    for e = 1:ne, m = 1:nk
        factor = im * energy_amplitude[e] * momentum_amplitude[m]
        lesser[e, m, :, :] .= factor .* (phased_v * phased_v')
        original_lesser[e, m, :, :] .= factor .* (v * v')
    end
    density = spatial_energy_density(problem, lesser)
    @test density.axis_order == (:energy, :position)
    @test size(density.n_per_eV_m3) == (ne, nz)
    @test eltype(density.n_per_eV_m3) == Float64
    @test density.z_m == problem.grids.x .* problem.scales.L₀_m
    @test density.energy_eV ==
          problem.grids.ε .* problem.scales.E₀_eV .+ ustrip(u"eV", problem.physical.E_ref)

    # Independent wavefunction oracle: rank-one density is |Σ_a χ_a v_a|²,
    # with no call to the matrix-projection implementation.
    wave_density = [abs2(sum(basis.χ[z, a]*v[a] for a = 1:nb)) for z = 1:nz]
    radial = dot(problem.grids.wᵏ, momentum_amplitude)
    prefactor = problem.physical.g_s / (2π*problem.scales.E₀_eV*problem.scales.L₀_m^3)
    oracle = prefactor .* (energy_amplitude * wave_density') .* radial
    @test density.n_per_eV_m3 ≈ oracle rtol=3e-14
    @test density.n_per_eV_m3 ≈
          spatial_energy_density(original, original_lesser).n_per_eV_m3 rtol=3e-14
    diagonal_only = [sum(abs2(basis.χ[z, a])*abs2(v[a]) for a = 1:nb) for z = 1:nz]
    @test norm(wave_density-diagonal_only) > 1e-3 * norm(wave_density)

    energy_weights_eV = problem.grids.wᴱ .* problem.scales.E₀_eV
    position_weights_m = problem.grids.wˣ .* problem.scales.L₀_m
    recovered_volume = density.n_per_eV_m3' * energy_weights_eV
    @test recovered_volume ≈ ustrip.(u"m^-3", electron_density(problem, lesser)) rtol=3e-14
    resolved_sheet = density.n_per_eV_m3 * position_weights_m
    oracle_sheet =
        problem.physical.g_s / (2π*problem.scales.E₀_eV*problem.scales.L₀_m^2) .*
        energy_amplitude .* radial
    @test resolved_sheet ≈ oracle_sheet rtol=3e-13
    integrated_sheet = dot(energy_weights_eV, resolved_sheet)
    @test integrated_sheet ≈
          real(tr(ustrip.(u"m^-2", sheet_density_matrix(problem, lesser)))) rtol=3e-13
    @test integrated_sheet ≈ sum(ustrip.(u"m^-2", state_populations(problem, lesser).Nₛ)) rtol=3e-13

    single_spin_physics = replace_fields(problem.physical, (g_s = 1, E_ref = 0.137u"eV"))
    single_spin = spatial_energy_density(
        replace_fields(problem, (physical = single_spin_physics,)),
        lesser,
    )
    @test single_spin.n_per_eV_m3 ≈ density.n_per_eV_m3 ./ problem.physical.g_s rtol=3e-14
    @test single_spin.energy_eV ≈ density.energy_eV .+ 0.137 rtol=3e-14

    # Historical spectral_maps has POSITION×ENERGY Unitful maps. Its N map
    # must use exactly the same spin, momentum and physical-unit factors.
    zeros4 = zeros(ComplexF64, size(lesser))
    green = GreenState(zeros4, lesser, zeros4, zeros4, ones(ne, nk), ones(ne, nk))
    historical = spectral_maps(problem, green)
    @test ustrip.(u"eV^-1*m^-3", historical.N) == permutedims(density.n_per_eV_m3)
    @test ustrip.(u"eV", historical.E) == density.energy_eV
    @test ustrip.(u"m", historical.z) == density.z_m

    @test_throws DimensionMismatch spatial_energy_density(
        problem,
        lesser[1:(end-1), :, :, :],
    )
    invalid = copy(lesser)
    invalid[1, 1, 1, 1] = ComplexF64(NaN, 0)
    @test_throws ArgumentError spatial_energy_density(problem, invalid)
    # Projection reports finite unphysical data honestly; no density clipping.
    @test spatial_energy_density(problem, -lesser).n_per_eV_m3 ≈ -density.n_per_eV_m3 rtol=3e-14
end

end # independent suite
