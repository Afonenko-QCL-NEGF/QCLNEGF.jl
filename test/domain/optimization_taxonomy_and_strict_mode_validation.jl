module Suite_T003
include("../support/common.jl")
include("../support/algorithm_modes.jl")

@testset "Optimization taxonomy and strict mode validation" begin
    exact = AlgorithmOptions(kernel_build = :direct)
    @test algorithm_impact(exact) == :physics_preserving
    @test algorithm_impact(
        AlgorithmOptions(
            contraction = :low_rank,
            low_rank_relative_tolerance = 1e-4,
            kernel_build = :direct,
        ),
    ) == :controlled_numerical
    @test algorithm_impact(
        AlgorithmOptions(
            contraction = :low_rank,
            low_rank_maximum_rank = 4,
            kernel_build = :direct,
        ),
    ) == :controlled_numerical
    @test algorithm_impact(
        AlgorithmOptions(localization = :none, kernel_build = :direct),
    ) == :controlled_numerical
    @test algorithm_impact(AlgorithmOptions(self_energy_structure = :diagonal)) ==
          :physical_model
    @test algorithm_manifest(exact)["hilbert"] == "fft"
    @test algorithm_manifest(exact)["anderson_history_depth"] == 5
    @test algorithm_manifest(exact)["anderson_damping"] == 1.0
    @test algorithm_manifest(exact)["anderson_regularization"] == 1e-12
    @test any(
        item -> item.id == :mode_space && item.status == :documented,
        optimization_catalog(),
    )
    @test_throws ArgumentError QCLNEGF._check_algorithm_options(
        AlgorithmOptions(hilbert = :cyclic),
    )
    @test_throws ArgumentError QCLNEGF._check_algorithm_options(
        AlgorithmOptions(low_rank_relative_tolerance = 1.0),
    )
    @test_throws ArgumentError QCLNEGF._check_algorithm_options(
        AlgorithmOptions(localization = :magic),
    )
end

end # independent suite
