module Suite_T085
include("../support/common.jl")
include("../support/production_parallelism.jl")

@testset "Production parallel scheduler and exact contractions" begin
    previous_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)
    try
        workers = min(4, Base.Threads.nthreads(:default))
        options = ProductionOptions(
            parallel_backend = :threads,
            worker_count = workers,
            verify_fft_roundoff = false,
        )
        trace = zeros(Int, 257)
        QCLNEGF._production_parallel_linear!(length(trace), options) do index
            trace[index] = Base.Threads.threadid()
        end
        @test all(>(0), trace)
        @test length(unique(trace)) == workers

        ranges = Vector{UnitRange{Int}}(undef, 8)
        QCLNEGF._production_worker_ranges!(49, 8) do worker, first_index, last_index
            ranges[worker] = first_index:last_index
        end
        @test all(!isempty, ranges)
        @test reduce(vcat, collect.(ranges)) == collect(1:49)
        @test maximum(length.(ranges)) - minimum(length.(ranges)) ≤ 1

        Random.seed!(0x50415241)
        NE, Nk, Nb = 19, 7, 3
        G = randn(ComplexF64, NE, Nk, Nb, Nb)
        K = randn(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
        wᵏ = rand(Nk)
        qᴷ = 0.071
        serial = production_static_contraction(
            K,
            G,
            wᵏ,
            qᴷ;
            energy_chunk = 4,
            parallel_backend = :blas,
        )
        threaded = production_static_contraction(
            K,
            G,
            wᵏ,
            qᴷ;
            energy_chunk = 4,
            parallel_backend = :threads,
            worker_count = workers,
        )
        # Every output column retains the same BLAS reduction; scheduling
        # changes only which independent energy chunk computes it.
        @test threaded == serial

        block = randn(ComplexF64, Nb, Nb, Nb, Nb)
        Kindependent = zeros(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
        for m = 1:Nk, mp = 1:Nk
            Kindependent[m, mp, :, :, :, :] .= block
        end
        serial_reduced = production_static_contraction(
            Kindependent,
            G,
            wᵏ,
            qᴷ;
            energy_chunk = 5,
            parallel_backend = :blas,
        )
        threaded_reduced = production_static_contraction(
            Kindependent,
            G,
            wᵏ,
            qᴷ;
            energy_chunk = 5,
            parallel_backend = :threads,
            worker_count = workers,
        )
        @test threaded_reduced == serial_reduced

        ε = collect(range(-0.8, 1.1; length = NE))
        h = randn(ComplexF64, Nk, Nb, Nb)
        for m = 1:Nk
            h[m, :, :] .= (h[m, :, :] + h[m, :, :]') / 2
        end
        Σᴿ = 0.03 .* randn(ComplexF64, NE, Nk, Nb, Nb)
        reference_GR, reference_condition, reference_scale =
            retarded_green(ε, h, Σᴿ; η = 2e-3, return_scale = true)
        production_GR, production_condition, production_scale =
            QCLNEGF._retarded_green_production(ε, h, Σᴿ; η = 2e-3, options)
        @test production_GR == reference_GR
        @test production_condition == reference_condition
        @test production_scale == reference_scale

        δ = 0.173
        plan = build_shift_plan(ε, δ)
        dense_shift = build_shift_matrix(ε, δ)
        T = randn(ComplexF64, Nb, Nb)
        embedding_reference =
            QCLNEGF._embedding_component(T, apply_energy_shift(dense_shift, G))
        embedding_threaded = QCLNEGF._embedding_component_production(T, G, plan, options)
        @test embedding_threaded ≈ embedding_reference rtol=2e-14 atol=2e-14
    finally
        BLAS.set_num_threads(previous_blas_threads)
    end
end

end # independent suite
