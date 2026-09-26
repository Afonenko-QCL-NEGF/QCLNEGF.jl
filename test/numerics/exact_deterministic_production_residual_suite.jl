module Suite_T087
include("../support/common.jl")
include("../support/production_residuals.jl")

@testset "Exact deterministic production residual suite" begin
    ε, h, green, total, candidate = _random_residual_fixture()
    workspace =
        ProductionResidualWorkspace(length(ε), size(h, 1), size(h, 2); chunk_size = 7)
    result = production_residual_suite(ε, h, green, total, candidate; workspace)

    reference = (
        r_D = QCLNEGF._dyson_residual(ε, h, total.Σᴿ, green.Gᴿ),
        r_A = QCLNEGF._spectral_residual(green, total),
        r_K = QCLNEGF._keldysh_residual(green, candidate),
        r_PSD = QCLNEGF._green_positivity(green, total),
        r_caus = QCLNEGF._causality_residual(total),
    )
    for name in keys(reference)
        @test getproperty(result, name) ≈ getproperty(reference, name) rtol=2e-14 atol=2e-14
    end

    # Repeated runs and different task partitioning have identical scalar
    # reductions; chunking changes scheduling, never the arithmetic order.
    @test production_residual_suite(ε, h, green, total, candidate; workspace) == result
    one_block_chunks =
        ProductionResidualWorkspace(length(ε), size(h, 1), size(h, 2); chunk_size = 1)
    whole_domain = ProductionResidualWorkspace(
        length(ε),
        size(h, 1),
        size(h, 2);
        chunk_size = length(ε) * size(h, 1),
    )
    @test production_residual_suite(
        ε,
        h,
        green,
        total,
        candidate;
        workspace = one_block_chunks,
    ) == result
    @test production_residual_suite(
        ε,
        h,
        green,
        total,
        candidate;
        workspace = whole_domain,
    ) == result
    single_worker = ProductionResidualWorkspace(
        length(ε),
        size(h, 1),
        size(h, 2);
        chunk_size = 7,
        worker_count = 1,
    )
    @test production_residual_suite(
        ε,
        h,
        green,
        total,
        candidate;
        workspace = single_worker,
    ) == result
    nested_result = fetch(
        Base.Threads.@spawn production_residual_suite(
            ε,
            h,
            green,
            total,
            candidate;
            workspace = single_worker,
        )
    )
    @test nested_result == result

    # All block products and LAPACK work arrays are preallocated.  The only
    # remaining allocations are bounded task bookkeeping, not `(N_E,N_k)`
    # matrix temporaries.
    production_residual_suite(ε, h, green, total, candidate; workspace)
    suite_bytes =
        @allocated production_residual_suite(ε, h, green, total, candidate; workspace)
    @test suite_bytes < 256_000

    # The convenience projections have the same exact result contract.
    @test QCLNEGF._dyson_residual_production(ε, h, green, total, candidate; workspace) ==
          result.r_D
    @test QCLNEGF._spectral_residual_production(
        ε,
        h,
        green,
        total,
        candidate;
        workspace,
    ) == result.r_A
    @test QCLNEGF._keldysh_residual_production(
        ε,
        h,
        green,
        total,
        candidate;
        workspace,
    ) == result.r_K
    @test QCLNEGF._green_positivity_production(
        ε,
        h,
        green,
        total,
        candidate;
        workspace,
    ) == result.r_PSD
    @test QCLNEGF._causality_residual_production(
        ε,
        h,
        green,
        total,
        candidate;
        workspace,
    ) == result.r_caus
end

end # independent suite
