module Suite_T091
include("../support/common.jl")

@testset "Threaded production seed equals the educational discretization" begin
    scattering = ScatteringOptions(
        LO = false,
        acoustic = false,
        impurity = false,
        IFR = false,
        alloy = false,
    )
    problem = build_problem(numerical = tutorial_numerics(), scattering = scattering)
    n = problem.numerical
    h = project_hamiltonians(problem, zeros(n.N_z))
    seed = QCLNEGF.effective_seed_parameters(problem)
    @test seed.requested_eV ≈ ustrip(u"eV", n.η_seed)
    @test seed.effective_eV == seed.requested_eV
    @test seed.grid_floor_eV ≈ ustrip(u"eV", (n.E_max-n.E_min)/(n.N_E-1))
    @test !seed.clamped
    @test seed.rule == :explicit_requested
    @test seed.working_eta_eV == 0.0

    reference, μ_reference = QCLNEGF._seed_green(problem, h)
    weight_threads = zeros(Int, n.N_E)
    lesser_threads = zeros(Int, n.N_E * n.N_k)
    production, μ_production = QCLNEGF._seed_green_production(
        problem,
        h;
        thread_trace = (weights = weight_threads, lesser = lesser_threads),
    )

    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    @test size(production.Gᴿ) == shape
    @test size(production.Gˡ) == shape
    @test size(production.Gᵍ) == shape
    @test size(production.A) == shape
    @test size(production.condition_number) == (n.N_E, n.N_k)
    @test size(production.dyson_scale) == (n.N_E, n.N_k)

    # Both paths use the same rescaled dense inversion; the production
    # primitive merely distributes independent (E,k) blocks over threads.
    @test production.Gᴿ == reference.Gᴿ
    @test production.A == reference.A
    @test production.condition_number == reference.condition_number
    @test production.dyson_scale == reference.dyson_scale
    @test μ_production ≈ μ_reference rtol=8eps(Float64) atol=8eps(Float64)
    @test production.Gˡ ≈ reference.Gˡ rtol=16eps(Float64) atol=16eps(Float64)
    @test production.Gᵍ ≈ reference.Gᵍ rtol=16eps(Float64) atol=16eps(Float64)

    target = QCLNEGF._scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    @test number_functional(production.Gˡ, problem.grids, problem.physical.g_s) ≈ target rtol=16eps(
        Float64,
    )

    # Every scheduled row/block must be visited.  With more than one Julia
    # worker and enough work, the static production loops must use the whole
    # default thread pool rather than silently reverting to one core.
    @test all(>(0), weight_threads)
    @test all(>(0), lesser_threads)
    expected_workers = min(Base.Threads.nthreads(:default), n.N_E)
    @test length(unique(weight_threads)) == expected_workers
    @test length(unique(lesser_threads)) ==
          min(Base.Threads.nthreads(:default), n.N_E * n.N_k)
end

end # independent suite
