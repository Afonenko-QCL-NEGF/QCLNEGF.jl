module SCBAPhysicalMarkerTests
include("../support/common.jl")
const BN = QCLNEGF
const NUM = BN.QCLNumerics

@testset "Sampled linewidth is a spectral modal energy-width diagnostic" begin
    energy = [-1.0, 0.0, 1.0]
    grids = (; ε = energy, wᴱ = [0.5, 1.0, 0.5], wᵏ = [1.0])
    rotation = ComplexF64[1 im; im 1]/sqrt(2)
    A = zeros(ComplexF64, 3, 1, 2, 2)
    sigmaR, sigmaL, sigmaG = (similar(A) for _ = 1:3)
    for e = 1:3
        A[e, 1, :, :] = rotation*Diagonal([1.0, 3.0])*rotation'
        gamma = rotation*Diagonal([0.1, 2.0])*rotation'
        sigmaR[e, 1, :, :] = -0.5im*gamma
        sigmaL[e, 1, :, :] = 0.4im*gamma
        sigmaG[e, 1, :, :] = -0.6im*gamma
    end
    total = SelfEnergyFamily(sigmaR, sigmaL, sigmaG)
    original = copy(A)
    policy = SCBAPhysicsMarkerPolicy(max_spectral_blocks = 3)
    marker = NUM._sampled_linewidth_markers(A, total, grids, policy)
    @test marker.status === :available
    @test marker.q10 ≈ 0.1
    @test marker.q50 ≈ 2.0
    @test marker.q90 ≈ 2.0
    @test marker.underresolved ≈ 0.25
    @test marker.blocks <= policy.max_spectral_blocks
    @test A == original
    # A change of energy units scales Γ and ΔE together; no angular-frequency factor.
    scaled = SelfEnergyFamily(7sigmaR, 7sigmaL, 7sigmaG)
    marker_scaled = NUM._sampled_linewidth_markers(
        A,
        scaled,
        (; ε = 7energy, wᴱ = 7grids.wᴱ, wᵏ = grids.wᵏ),
        policy,
    )
    @test marker_scaled.q50 ≈ marker.q50
    # Sub-roundoff negative modes are ignored in diagnostics only.
    for e = 1:3
        A[e, 1, :, :] = Diagonal(ComplexF64[1.0, -1e-18])
    end
    @test NUM._sampled_linewidth_markers(A, total, grids, policy).status === :available
    @test real(A[1, 1, 2, 2]) == -1e-18
    A[:, 1, 2, 2] .= -0.5
    @test NUM._sampled_linewidth_markers(A, total, grids, policy).status ===
          :indefinite_spectrum
end

@testset "Raw and normalized matrix FDT differ and retain the spectral chemical potential" begin
    physical = reference_parameters(V_period = 0u"mV", Tᴸ = 70u"K")
    scattering = ScatteringOptions(
        LO = false,
        acoustic = false,
        impurity = false,
        IFR = false,
        alloy = false,
    )
    problem = build_problem(; physical, numerical = tutorial_numerics(), scattering)
    h = project_hamiltonians(problem, zeros(problem.numerical.N_z))
    n = problem.numerical
    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    sigmaR, sigmaL, sigmaG = (zeros(ComplexF64, shape) for _ = 1:3)
    gamma = 0.02
    for e = 1:n.N_E, m = 1:n.N_k, a = 1:n.N_b
        sigmaR[e, m, a, a] = -0.5im*gamma
    end
    GR, _, _ = BN._retarded_green_production(problem.grids.ε, h, sigmaR)
    A = spectral_function(GR)
    weights = [
        sum(problem.grids.wᵏ[m]*real(tr(A[e, m, :, :])) for m = 1:n.N_k) *
        problem.grids.wᴱ[e] *
        physical.g_s/(2π) for e = 1:n.N_E
    ]
    target = BN._scaled_physics(physical, problem.scales).Nᴰ²ᴰ
    kBT = ustrip(u"eV", BN.CODATA.kᴮₑᵥ*physical.Tᴸ)/problem.scales.E₀_eV
    fit = NUM._marker_spectral_chemical_potential(problem.grids.ε, weights, target, kBT)
    @test fit.status === :available
    @test sum(weights[e]*NUM._fermi(problem.grids.ε[e], fit.mu, kBT) for e = 1:n.N_E) ≈
          target rtol=1e-12
    @test NUM._marker_spectral_chemical_potential(
        problem.grids.ε,
        weights,
        2sum(weights),
        kBT,
    ).status === :spectral_number_not_bracketed
    for deficit in (1.0, 0.99)
        for e = 1:n.N_E, m = 1:n.N_k, a = 1:n.N_b
            f = deficit*NUM._fermi(problem.grids.ε[e], fit.mu, kBT)
            sigmaL[e, m, a, a] = im*f*gamma
            sigmaG[e, m, a, a] = -im*(1-f)*gamma
        end
        embedding = SelfEnergyFamily(sigmaR, sigmaL, sigmaG)
        green, lambda, total = BN._current_green_production(
            problem,
            h,
            Dict{Symbol,SelfEnergyFamily}(),
            embedding,
            ProductionOptions(),
        )
        marker = NUM._scba_physical_markers(
            problem,
            green,
            total,
            -3.0,
            7,
            SCBAPhysicsMarkerPolicy(max_spectral_blocks = 32),
        )
        source_marker = NUM._scba_physical_markers(
            problem,
            green,
            total,
            -3.0,
            7,
            SCBAPhysicsMarkerPolicy(max_spectral_blocks = 32);
            seed_mu_scaled = fit.mu,
            seed_number_ratio = 1.0,
        )
        @test source_marker.seed_mu_eV ≈ fit.mu*problem.scales.E₀_eV
        @test source_marker.seed_number_ratio == 1.0
        @test source_marker.fdt_raw_seed_mu ≈ marker.fdt_raw atol=1e-13
        @test source_marker.fdt_normalized_seed_mu ≈ marker.fdt_normalized atol=1e-13
        shifted_marker = NUM._scba_physical_markers(
            problem,
            green,
            total,
            -3.0,
            7,
            SCBAPhysicsMarkerPolicy(max_spectral_blocks = 32);
            seed_mu_scaled = fit.mu+2kBT,
            seed_number_ratio = 1.0,
        )
        @test shifted_marker.fdt_raw == source_marker.fdt_raw
        @test shifted_marker.fdt_raw_seed_mu > source_marker.fdt_raw_seed_mu
        @test isnan(marker.seed_mu_eV)
        @test isnan(marker.fdt_raw_seed_mu)
        @test marker.equilibrium_applicable
        @test marker.equilibrium_status === :available
        @test marker.equilibrium_abs_current_A_m2 == 3.0
        @test marker.equilibrium_mu_eV ≈ fit.mu*problem.scales.E₀_eV
        @test marker.represented_capacity ≈ sum(weights) rtol=1e-12
        @test marker.raw_hole_charge ≈ sum(weights)-deficit*target rtol=1e-12
        @test marker.occupied_fraction_a ≈ 1.0 atol=2e-12
        @test marker.empty_fraction_c ≈ (1-deficit)*target/marker.raw_hole_charge atol=1e-14
        if deficit == 1.0
            @test marker.fdt_raw < 1e-12
            @test marker.fdt_normalized < 1e-12
            @test marker.relative_correction_Gn < 1e-11
        else
            @test marker.fdt_raw > 0
            @test marker.fdt_normalized > 0
            @test marker.relative_correction_Gn > 0
            @test marker.relative_correction_Gp > 0
            # Independent matrix reconstruction, including off-diagonal entries.
            raw_error = normalized_error = scale = 0.0
            for e = 1:n.N_E, m = 1:n.N_k
                f = NUM._fermi(problem.grids.ε[e], fit.mu, kBT)
                rawGn = -im*green.Gᴿ[e, m, :, :]*total.Σˡ[e, m, :, :]*green.Gᴿ[e, m, :, :]'
                rawGp = im*green.Gᴿ[e, m, :, :]*total.Σᵍ[e, m, :, :]*green.Gᴿ[e, m, :, :]'
                w = problem.grids.wᴱ[e]*problem.grids.wᵏ[m]
                raw_error +=
                    w*(norm(rawGn-f*A[e, m, :, :])^2+norm(rawGp-(1-f)*A[e, m, :, :])^2)
                normalized_error +=
                    w*(
                        norm(-im*green.Gˡ[e, m, :, :]-f*A[e, m, :, :])^2 +
                        norm(im*green.Gᵍ[e, m, :, :]-(1-f)*A[e, m, :, :])^2
                    )
                scale += w*norm(A[e, m, :, :])^2
            end
            @test marker.fdt_raw ≈ sqrt(raw_error/scale) rtol=2e-11
            @test marker.fdt_normalized ≈ sqrt(normalized_error/scale) rtol=2e-11
        end
    end
    biased = build_problem(numerical = tutorial_numerics(), scattering = scattering)
    @test NUM._marker_equilibrium_status(biased) === :finite_bias
end

@testset "Marker cadence never loses a terminal or checkpoint observation" begin
    problem = build_problem(numerical = tutorial_numerics())
    events = BN.SolverEvent[]
    observed = SCBAIteration[]
    checkpoints = SCBAResult[]
    options = SolverOptions(max_scba = 4, max_poisson = 1)
    production = ProductionOptions(
        phase_timing = true,
        physics_markers = SCBAPhysicsMarkerPolicy(cadence = 10, max_spectral_blocks = 8),
        checkpoint_every_scba = 2,
        event_sink = event->push!(events, event),
    )
    result = solve_scba_production(
        problem,
        zeros(problem.numerical.N_z);
        options,
        production_options = production,
        history_observer = row->push!(observed, row),
        iteration_callback = row->push!(checkpoints, row),
    )
    @test length(observed) == length(result.history)
    @test first(result.history).physical_markers !== nothing
    @test isfinite(first(result.history).physical_markers.seed_mu_eV)
    @test first(result.history).physical_markers.seed_number_ratio ≈ 1.0 rtol=1e-10
    @test length(result.history) == 4
    @test result.history[3].physical_markers === nothing
    @test last(result.history).physical_markers !== nothing
    @test last(result.history).physical_markers.measured_iteration == last(result.history).ν
    @test last(result.history).physical_markers.fresh_map_status === :available
    @test length(last(result.history).physical_markers.collisions) ==
          2length(problem.kernels.enabled)
    @test all(
        row.state_kind in (:mixed, :fresh) for
        row in last(result.history).physical_markers.collisions
    )
    @test all(
        last(checkpoint.history).physical_markers !== nothing for checkpoint in checkpoints
    )
    timed = [
        Dict(metric.name=>metric.value for metric in event.metrics) for
        event in events if any(metric.name===:t_candidate for metric in event.metrics)
    ]
    @test !isempty(timed)
    for sample in timed
        @test sample[:t_candidate] >=
              sample[:t_candidate_scattering]+sample[:t_candidate_embedding]+sample[:t_candidate_assembly]
        @test sample[:t_physics_markers] >= 0
    end
    # Native v4 fixtures explicitly declare unavailable marker/witness slots.
    row = SCBAIteration(
        1,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        NaN,
        NaN,
        NaN,
        NaN,
        nothing,
        nothing,
    )
    @test row.physical_markers === nothing
    @test_throws MethodError SCBAIteration(1, zeros(12)...)
end

@testset "Fresh channel components retain raw scale and same-state collision balances" begin
    shape = (3, 1, 2, 2)
    old = SelfEnergyFamily(fill(1.0im, shape), fill(2.0im, shape), fill(-3.0im, shape))
    fresh = SelfEnergyFamily(fill(2.0im, shape), fill(3.0im, shape), fill(-5.0im, shape))
    rows = BN.SCBAChannelMarker[]
    NUM._channel_component_markers!(rows, :LO, old, fresh)
    for (row, field) in zip(rows, (:Σᴿ, :Σˡ, :Σᵍ))
        @test row.residual_absolute ≈ norm(getfield(fresh, field)-getfield(old, field))
        @test row.residual_scale ≈ norm(getfield(fresh, field))
        @test row.residual_relative ≈ row.residual_absolute/(row.residual_scale+1e-14)
    end
    green =
        (; Gᴿ = zeros(ComplexF64, shape), Gˡ = fill(0.2im, shape), Gᵍ = fill(-0.8im, shape))
    problem = (;
        grids = (; ε = [-1.0, 0.0, 1.0], wᴱ = [0.5, 1.0, 0.5], wᵏ = [1.0]),
        physical = (; g_s = 2),
    )
    observed = NUM._channel_collision_marker(problem, green, :LO, :fresh, fresh)
    expected =
        sum(
            problem.grids.wᴱ[e]*tr(
                fresh.Σˡ[e, 1, :, :]*green.Gᵍ[e, 1, :, :]-fresh.Σᵍ[e, 1, :, :]*green.Gˡ[
                    e,
                    1,
                    :,
                    :,
                ],
            ) for e = 1:3
        )/π
    @test observed.particle_signed ≈ real(expected)
    @test observed.particle_imaginary ≈ imag(expected)
    @test observed.energy_signed ≈ 0 atol=1e-14
    @test observed.particle_absolute >= abs(observed.particle_signed)
    @test observed.energy_absolute > 0
    mixed = NUM._channel_collision_marker(problem, green, :LO, :mixed, old)
    @test mixed.particle_signed != observed.particle_signed
end

@testset "Shift boundaries expose represented weight and never invent exterior flux" begin
    g = (;
        Gᴿ = zeros(ComplexF64, 3, 1, 1, 1),
        Gˡ = fill(1im, 3, 1, 1, 1),
        A = fill(2.0, 3, 1, 1, 1),
    )
    p = (; grids = (; ε = [-1.0, 0.0, 1.0], wᴱ = [0.5, 1.0, 0.5], wᵏ = [1.0]))
    @test NUM._shift_boundary_markers(p, g, 0.0) == (0.0, 0.0)
    @test NUM._shift_boundary_markers(p, g, 0.5) == (0.5, 0.5)
    @test NUM._shift_boundary_markers(p, g, 3.0) == (1.0, 1.0)
end


@testset "Threshold summaries separate transient crossing from sustained passage" begin
    row(i, r) = SCBAIteration(
        i,
        0.0,
        0.0,
        r,
        r,
        r,
        1.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        100.0+i,
        1.0,
        1.0,
        1.0,
        0.0,
        nothing,
        nothing,
    )
    history = [row(i, r) for (i, r) in enumerate([1e-2, 1e-5, 1e-2, 1e-5, 1e-5, 1e-5])]
    summaries = NUM._scba_threshold_crossings(history; required_consecutive = 3)
    @test length(summaries) == 20
    crossing = only(x for x in summaries if x.metric === :joint && x.threshold == 1e-4)
    @test crossing.first_iteration == 2
    @test crossing.sustained_start_iteration == 4
    @test crossing.sustained_end_iteration == 6
    @test crossing.first_current_A_m2 == 102.0
    @test crossing.sustained_current_A_m2 == 106.0
    missing = only(x for x in summaries if x.metric === :joint && x.threshold == 1e-8)
    @test missing.first_iteration == 0
    @test isnan(missing.first_current_A_m2)
    gaps = NUM._scba_threshold_crossings([row(1, 1e-5), row(3, 1e-5), row(4, 1e-5)])
    @test all(x.sustained_end_iteration == 0 for x in gaps)
    @test_throws ArgumentError NUM._scba_threshold_crossings([row(1, 1e-5), row(1, 1e-5)])
    @test all(
        x.first_iteration == 0 for x in NUM._scba_threshold_crossings(SCBAIteration[])
    )
end

end
