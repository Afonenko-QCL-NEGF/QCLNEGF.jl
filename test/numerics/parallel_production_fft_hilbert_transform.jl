module Suite_T082
include("../support/common.jl")
include("../support/production_fft_parallel.jl")

@testset "Parallel production FFT Hilbert transform" begin
    Random.seed!(0xB05C0FF7)
    NE, Nk, Nb = 41, 3, 2
    ε = collect(range(-1.1, 1.4; length = NE))
    Δε = ε[2] - ε[1]
    wᴱ = fill(Δε, NE)
    wᴱ[[1, end]] ./= 2
    Γ = randn(ComplexF64, NE, Nk, Nb, Nb)
    columns = Nk * Nb^2
    requested_workers = min(4, Base.Threads.nthreads(:default))

    planner_threads_before = FFTW.get_num_threads()
    plan = ProductionFFTHilbertPlan(
        ε,
        columns;
        column_chunk = 3,
        worker_count = requested_workers,
        verify_roundoff = true,
    )
    @test FFTW.get_num_threads() == planner_threads_before

    ranges = getfield.(plan.primary.workers, :columns)
    @test reduce(vcat, collect.(ranges)) == collect(1:columns)
    @test all(
        isempty(intersect(ranges[i], ranges[j])) for
        i in eachindex(ranges), j in eachindex(ranges) if i < j
    )
    @test length(ranges) == min(requested_workers, columns)
    expected_buffer_columns = min(plan.column_chunk, maximum(length, ranges))
    @test all(
        size(worker.buffer) == (plan.primary.transform_length, expected_buffer_columns) for
        worker in plan.primary.workers
    )

    direct = direct_hilbert_transform(Γ, ε, wᴱ)
    serial = fft_hilbert_transform(Γ, ε, wᴱ; column_chunk = 3)
    parallel, rround =
        production_fft_hilbert_transform(Γ, ε, wᴱ, plan; return_roundoff = true)
    @test parallel ≈ serial rtol=2e-14 atol=2e-14
    @test parallel ≈ direct rtol=2e-13 atol=2e-13
    @test isfinite(rround)
    @test rround < 1e-12

    # Reusing the same plans must neither retain input data in inactive buffer
    # columns nor change the result.
    Γsecond = randn(ComplexF64, NE, Nk, Nb, Nb)
    reused = production_fft_hilbert_transform(Γsecond, ε, wᴱ, plan)
    serial_second = fft_hilbert_transform(Γsecond, ε, wᴱ; column_chunk = 3)
    @test reused ≈ serial_second rtol=2e-14 atol=2e-14

    one_shot = production_fft_hilbert_transform(
        Γ,
        ε,
        wᴱ;
        column_chunk = 5,
        worker_count = requested_workers,
    )
    @test one_shot ≈ serial rtol=2e-14 atol=2e-14

    C = sizeof(ComplexF64)
    expected_primary =
        plan.primary.transform_length *
        C *
        (1 + length(plan.primary.workers) * expected_buffer_columns)
    verification_ranges = getfield.(plan.verification.workers, :columns)
    expected_verification_columns =
        min(plan.column_chunk, maximum(length, verification_ranges))
    expected_verification =
        plan.verification.transform_length *
        C *
        (1 + length(plan.verification.workers) * expected_verification_columns)
    @test production_fft_workspace_bytes(plan) == expected_primary
    @test production_fft_workspace_bytes(plan; verification = true) == expected_verification

    primary_only = ProductionFFTHilbertPlan(
        ε,
        columns;
        column_chunk = 3,
        worker_count = requested_workers,
    )
    @test_throws ArgumentError production_fft_hilbert_transform(
        Γ,
        ε,
        wᴱ,
        primary_only;
        return_roundoff = true,
    )
    @test_throws DimensionMismatch production_fft_hilbert_transform(
        Γ[:, 1:2, :, :],
        ε,
        wᴱ,
        plan,
    )
    @test_throws ArgumentError ProductionFFTHilbertPlan(ε, columns; worker_count = 0)

    # The admission estimator must account for the exact fixed-range buffer
    # geometry used above: one kernel vector plus one bounded buffer per active
    # worker, and an additional doubled-length workspace only for verification.
    numerical = tutorial_numerics()
    estimated_columns = numerical.N_k * numerical.N_b^2
    estimated_workers = min(requested_workers, estimated_columns)
    estimated_buffer_columns = min(3, cld(estimated_columns, estimated_workers))
    estimated_length = nextpow(2, 2 * numerical.N_E - 1)
    expected_estimated_fft =
        3 *
        estimated_length *
        sizeof(ComplexF64) *
        (1 + estimated_workers * estimated_buffer_columns)
    memory_options = ProductionOptions(
        hilbert_columns = 3,
        worker_count = requested_workers,
        verify_fft_roundoff = true,
    )
    memory_estimate = estimate_production_memory(
        numerical,
        1;
        dense_mechanism_count = 1,
        lo_kernel = :dense,
        options = memory_options,
    )
    @test memory_estimate.breakdown[:fft_workspace] == expected_estimated_fft

    primary_memory_options = ProductionOptions(
        hilbert_columns = 3,
        parallel_backend = :blas,
        worker_count = requested_workers,
        verify_fft_roundoff = false,
    )
    primary_memory_estimate = estimate_production_memory(
        numerical,
        1;
        dense_mechanism_count = 1,
        lo_kernel = :dense,
        options = primary_memory_options,
    )
    @test primary_memory_estimate.breakdown[:fft_workspace] == expected_estimated_fft ÷ 3
end

end # independent suite
