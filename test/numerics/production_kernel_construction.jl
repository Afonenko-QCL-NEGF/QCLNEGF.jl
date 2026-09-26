module Suite_T084
include("../support/common.jl")

@testset "Production kernel construction" begin
    p = reference_parameters()
    n = tutorial_numerics()
    s = ScaleSystem()
    grids = build_grids(p, n, s)
    profiles = build_profiles(p, n, s, grids)
    basis = build_localized_basis(p, n, s, grids, profiles)
    scattering = default_scattering()

    direct = build_kernels(p, n, scattering, s, grids, profiles, basis)
    exact, exact_diagnostic =
        QCLNEGF.build_kernels_exact_parallel(p, n, scattering, s, grids, profiles, basis)
    production, diagnostic = QCLNEGF.build_kernels_production(
        p,
        n,
        scattering,
        s,
        grids,
        profiles,
        basis;
        options = QCLNEGF.ProductionKernelOptions(
            initial_lookup_nodes = 33,
            maximum_lookup_nodes = 257,
            relative_tolerance = 5e-4,
            angular_validation_pairs = 3,
            strict = true,
        ),
    )

    @test exact.enabled == direct.enabled
    @test exact_diagnostic.all_accepted
    @test exact_diagnostic.worst_measured_relative_residual == 0.0
    @test exact.Fᴸᴼ ≈ direct.Fᴸᴼ rtol=8eps(Float64) atol=8eps(Float64)
    for mechanism in direct.enabled
        raw_direct = direct.qᴷ[mechanism] .* direct.K[mechanism]
        raw_exact = exact.qᴷ[mechanism] .* exact.K[mechanism]
        residual = norm(raw_exact - raw_direct) / max(norm(raw_direct), eps(Float64))
        # The exact production builder changes only independent block order
        # and dense basis contractions. It has no lookup/interpolation path.
        @test residual < 2e-12
        exact_record = exact_diagnostic.mechanisms[mechanism]
        @test exact_record.accepted
        @test exact_record.lookup_nodes == 0
        @test exact_record.direct_q_evaluations == 0
        @test startswith(String(exact_record.construction), "exact_")
    end

    @test production.enabled == direct.enabled
    @test diagnostic.all_accepted
    @test diagnostic.worst_measured_relative_residual <= 5e-4
    @test production.Fᴸᴼ ≈ direct.Fᴸᴼ rtol=8eps(Float64) atol=8eps(Float64)

    for mechanism in direct.enabled
        raw_direct = direct.qᴷ[mechanism] .* direct.K[mechanism]
        raw_production = production.qᴷ[mechanism] .* production.K[mechanism]
        residual = norm(raw_production - raw_direct) / norm(raw_direct)
        if mechanism in (:acoustic, :IFR)
            @test residual < 2e-13
            @test diagnostic.mechanisms[mechanism].construction !=
                  :piecewise_linear_q_lookup
        else
            # This is an independent full-tensor comparison.  It is kept
            # somewhat looser than the sampled entrywise table criterion
            # because it measures a different norm after angular averaging.
            @test residual < 2e-3
            d = diagnostic.mechanisms[mechanism]
            @test d.construction == :piecewise_linear_q_lookup
            @test d.direct_q_evaluations > d.lookup_nodes
            @test length(d.node_history) == length(d.residual_history)
            @test d.accepted
        end
        rH, rPSD = QCLNEGF._kernel_residual(production.K[mechanism])
        @test rH < 2e-10
        @test rPSD < 2e-10
    end

    @test_throws ArgumentError QCLNEGF.ProductionKernelOptions(
        initial_lookup_nodes = 2,
    ) |> QCLNEGF._check_production_kernel_options

    raw_fixture = Dict(:acoustic => fill(2.0 + 0.0im, 1, 1, 1, 1, 1, 1))
    raw_storage = raw_fixture[:acoustic]
    normalized_fixture =
        QCLNEGF._normalize_production_kernels(raw_fixture, zeros(ComplexF64, 1, 1, 1))
    @test normalized_fixture.K[:acoustic] === raw_storage
    @test only(normalized_fixture.K[:acoustic]) == 1.0 + 0.0im
    @test normalized_fixture.qᴷ[:acoustic] == 2.0
end

end # independent suite
