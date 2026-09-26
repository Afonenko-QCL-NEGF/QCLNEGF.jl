module Suite_T072
include("../support/common.jl")

@testset "Physics-changing modes satisfy their declared contracts" begin
    problem = build_problem(
        numerical = tutorial_numerics(),
        scattering = ScatteringOptions(
            LO = false,
            acoustic = true,
            impurity = false,
            IFR = false,
            alloy = false,
        ),
    )
    h = project_hamiltonians(problem, zeros(problem.numerical.N_z))
    options = ProductionOptions(
        progress_every_scba = 0,
        progress_every_outer = 0,
        verify_fft_roundoff = false,
        algorithms = AlgorithmOptions(
            kernel_build = :direct,
            retarded_real_part = :drop,
            self_energy_structure = :diagonal,
            transverse_momentum = :averaged,
        ),
    )
    green, _ = QCLNEGF._seed_green_production(problem, h, options)
    cache = build_production_cache(problem; options)
    candidate =
        QCLNEGF._scattering_candidate_production(problem, green, cache, options)[:acoustic]
    @test algorithm_impact(options.algorithms) == :physical_model
    @test all(iszero, candidate.Σᴿ[:, :, 1, 2])
    @test all(iszero, candidate.Σˡ[:, :, 2, 1])
    for component in (:Σᴿ, :Σˡ, :Σᵍ), m = 2:problem.numerical.N_k
        @test getfield(candidate, component)[:, m, :, :] ==
              getfield(candidate, component)[:, 1, :, :]
    end
    Γ = im .* (candidate.Σᵍ .- candidate.Σˡ)
    @test candidate.Σᴿ ≈ -0.5im .* Γ rtol=8eps(Float64) atol=8eps(Float64)
end

end # independent suite
