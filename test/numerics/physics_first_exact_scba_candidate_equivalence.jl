module Suite_T071
include("../support/common.jl")

@testset "Physics-first exact SCBA candidate equivalence" begin
    problem = build_problem(numerical = tutorial_numerics())
    hartree = zeros(problem.numerical.N_z)
    h = project_hamiltonians(problem, hartree)
    seed_options = ProductionOptions(
        progress_every_scba = 0,
        progress_every_outer = 0,
        verify_fft_roundoff = false,
    )
    green, _ = QCLNEGF._seed_green_production(problem, h, seed_options)

    literal_options = ProductionOptions(
        progress_every_scba = 0,
        progress_every_outer = 0,
        algorithms = AlgorithmOptions(
            solver_backend = :production,
            energy_shift = :dense,
            hilbert = :direct,
            contraction = :literal,
            kernel_build = :direct,
            mixing = :linear,
        ),
    )
    optimized_options = ProductionOptions(
        progress_every_scba = 0,
        progress_every_outer = 0,
        algorithms = AlgorithmOptions(
            solver_backend = :production,
            energy_shift = :sparse_plan,
            hilbert = :fft,
            contraction = :dense_blas,
            kernel_build = :direct,
            mixing = :linear,
        ),
    )
    literal_cache = build_production_cache(problem; options = literal_options)
    optimized_cache = build_production_cache(problem; options = optimized_options)
    literal, _ = QCLNEGF._scattering_candidate_production(
        problem,
        green,
        literal_cache,
        literal_options;
        return_roundoff = true,
    )
    optimized, roundoff = QCLNEGF._scattering_candidate_production(
        problem,
        green,
        optimized_cache,
        optimized_options;
        return_roundoff = true,
    )
    @test roundoff < 1e-11
    for mechanism in problem.kernels.enabled
        for component in (:Σᴿ, :Σˡ, :Σᵍ)
            reference = getfield(literal[mechanism], component)
            candidate = getfield(optimized[mechanism], component)
            relative = norm(candidate - reference) / max(norm(reference), 1e-14)
            @test relative < 2e-11
        end
    end

    dense_embedding = QCLNEGF._embedding_self_energy_selected(
        problem,
        green,
        literal_cache,
        literal_options,
    )
    sparse_embedding = QCLNEGF._embedding_self_energy_selected(
        problem,
        green,
        optimized_cache,
        optimized_options,
    )
    for family_index = 1:3, component in (:Σᴿ, :Σˡ, :Σᵍ)
        @test getfield(sparse_embedding[family_index], component) ≈
              getfield(dense_embedding[family_index], component) rtol=3e-14 atol=3e-14
    end
end

end # independent suite
