module Suite_T111
include("../support/common.jl")

@testset "full-window basis preserves the spectral subspace" begin
    p, s = reference_parameters(), ScaleSystem()
    n = NumericalParameters(
        N_z = 18,
        N_b = 3,
        P_basis = 1,
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
    grids = build_grids(p, n, s)
    profiles = build_profiles(p, n, s, grids)
    for method in (:pzp_tails, :bloch_wannier, :wannier_stark)
        reference = QCLNEGF.build_multiperiod_reference_basis(
            p,
            n,
            s,
            grids,
            profiles;
            localization = method,
        )
        @test reference.periods == 3
        @test size(reference.transform) == (54, 9)
        @test reference.overlap ≈ I atol=2e-12
        @test reference.subspace_residual < 1e-12
        @test eigvals(Hermitian(reference.hamiltonian)) ≈ reference.reference_energies rtol=2e-12 atol=2e-12
        @test reference.hamiltonian ≈ QCLNEGF.project_reference_operator(
            reference,
            reference.real_space_hamiltonian,
        ) atol=2e-13
        @test all(reference.central_tail_weights .>= -1e-13)
        @test reference.coupling_truncation == :none
        # Complete coherent projection of a spatial density operator preserves
        # the expectation value in the original window, including all tails.
        state = ComplexF64.(1:9)
        state ./= norm(state)
        density_operator = Diagonal(cos.(reference.coordinates))
        ψ = reference.transform * state
        @test dot(ψ, density_operator*ψ) ≈ dot(
            state,
            QCLNEGF.project_reference_operator(reference, density_operator)*state,
        ) atol=1e-13
    end
    @test_throws ArgumentError build_basis(
        p,
        n,
        s,
        grids,
        profiles;
        localization = :real_space,
    )
    all_nodes = Dict(name => getfield(n, name) for name in fieldnames(NumericalParameters))
    all_nodes[:N_b] = n.N_z
    full_numerical = NumericalParameters(; all_nodes...)
    exact_basis =
        build_basis(p, full_numerical, s, grids, profiles; localization = :real_space)
    @test exact_basis.Φ ≈ I
    dispersion = QCLNEGF.basis_dispersion_diagnostic(exact_basis, profiles, grids, s)
    @test dispersion.maximum_error_eV < 1e-11

    @test_throws ArgumentError build_basis(
        p,
        n,
        s,
        grids,
        profiles;
        localization = :legacy_pzp,
    )
    canonical_basis = build_basis(p, n, s, grids, profiles; localization = :pzp)

    # The independent literal vertex contraction and optimized BLAS evaluator
    # are integrated by the same continuum angular rule; they must agree.
    literal = QCLNEGF.impurity_kernel(
        p,
        s,
        grids,
        profiles,
        canonical_basis;
        angular = :adaptive,
        angular_tolerance = 1e-5,
    )
    scattering = ScatteringOptions(
        LO = false,
        acoustic = false,
        impurity = true,
        IFR = false,
        alloy = false,
    )
    optimized, diagnostics = QCLNEGF.build_kernels_exact_parallel(
        p,
        n,
        scattering,
        s,
        grids,
        profiles,
        canonical_basis;
        impurity_angular = :adaptive,
        angular_tolerance = 1e-5,
    )
    @test optimized.qᴷ[:impurity] .* optimized.K[:impurity] ≈ literal rtol=2e-11
    @test diagnostics.mechanisms[:impurity].angular_quadrature_checked
    @test diagnostics.mechanisms[:impurity].angular_quadrature_relative_error !== nothing

    # Configured cavity adapter must use the same declared energy remap for
    # all three injection components, including fractional field shifts.
    empty_scattering = ScatteringOptions(
        LO = false,
        acoustic = false,
        impurity = false,
        IFR = false,
        alloy = false,
    )
    problem = build_problem(
        physical = p,
        numerical = n,
        scattering = empty_scattering,
        scales = s,
        validate_static = false,
    )
    h = project_hamiltonians(problem, zeros(n.N_z))
    R = zeros(ComplexF64, n.N_E, n.N_k, n.N_b, n.N_b)
    L, G = copy(R), copy(R)
    for e = 1:n.N_E, m = 1:n.N_k, a = 1:n.N_b
        gamma, occupation = 0.09+0.001e, 0.2+0.005e
        R[e, m, a, a] = -0.5im*gamma
        L[e, m, a, a] = im*gamma*occupation
        G[e, m, a, a] = -im*gamma*(1-occupation)
    end
    family = SelfEnergyFamily(R, L, G)
    _, plus, minus = QCLNEGF.cavity_embedding_self_energy(
        problem,
        family,
        h;
        periods = 1,
        energy_shift = :conservative_pair,
    )
    e, m = 16, 2
    drop = QCLNEGF._scaled_physics(p, s).Eᵖ
    shift = QCLNEGF.build_conservative_shift_pair(grids.ε, grids.wᴱ, drop)
    shifts = (shift.minus, Matrix{Float64}(I, n.N_E, n.N_E), shift.plus)
    blocks(field) =
        [Matrix(view(QCLNEGF.apply_energy_shift(W, field), e, m, :, :)) for W in shifts]
    onsite =
        [Matrix(view(h, m, :, :))-j*drop*Matrix{ComplexF64}(I, n.N_b, n.N_b) for j = -1:1]
    reference = QCLNEGF.finite_chain_green(
        grids.ε[e],
        onsite,
        [problem.basis.T₊, problem.basis.T₊];
        sigma_retarded = blocks(R),
        sigma_lesser = blocks(L),
        sigma_greater = blocks(G),
    )
    for (field, component) in ((:Σᴿ, :retarded), (:Σˡ, :lesser), (:Σᵍ, :greater))
        @test Matrix(view(getproperty(plus, field), e, m, :, :)) ≈
              getproperty(reference.plus, component) rtol=2e-12
        @test Matrix(view(getproperty(minus, field), e, m, :, :)) ≈
              getproperty(reference.minus, component) rtol=2e-12
    end
end

end # independent suite
