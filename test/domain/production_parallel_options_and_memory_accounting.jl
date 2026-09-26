module Suite_T086
include("../support/common.jl")
include("../support/production_parallelism.jl")

@testset "Production parallel options and memory accounting" begin
    @test_throws ArgumentError QCLNEGF._check_production_options(
        ProductionOptions(parallel_backend = :invalid),
    )
    @test_throws ArgumentError QCLNEGF._check_production_options(
        ProductionOptions(worker_count = -1),
    )
    @test_throws ArgumentError QCLNEGF._check_production_options(
        ProductionOptions(residual_chunk = 0),
    )
    @test_throws ArgumentError QCLNEGF._production_worker_count(:threads, -1, 10)
    @test_throws ArgumentError QCLNEGF._production_worker_count(:threads, 0, -1)
    @test_throws ArgumentError QCLNEGF._production_worker_count(:invalid, 0, 0)

    n = reference_production_numerics()
    options = ProductionOptions(
        parallel_backend = :threads,
        worker_count = min(4, Base.Threads.nthreads(:default)),
        energy_chunk = 64,
        hilbert_columns = 32,
        residual_chunk = 256,
    )
    estimate = estimate_production_memory(
        n,
        4;
        dense_mechanism_count = 3,
        lo_kernel = :dense,
        options,
    )
    Ncompound = n.N_k * n.N_b^2
    expected_work = 8.0 * n.N_E * (8 * Ncompound^2 + 2 * (n.N_k * n.N_b^2 + n.N_b^4))
    @test estimate.kernel_flops_per_candidate == expected_work
    @test estimate.breakdown[:fft_workspace] > 0
    @test estimate.breakdown[:residual_workspace] > 0
    @test estimate.breakdown[:scalar_reduction_workspace] > 0
    @test estimate.peak_bytes < options.memory_budget_bytes
end

end # independent suite
