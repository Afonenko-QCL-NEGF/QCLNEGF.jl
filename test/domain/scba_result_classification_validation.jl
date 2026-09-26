module Suite_T060
include("../support/common.jl")
include("../support/numerical_convergence.jl")

@testset "SCBA result classification validation" begin
    @test QCLNEGF._check_scba_result_classification(
        true,
        :converged,
        :strictly_converged,
    ) == :strictly_converged
    @test QCLNEGF._check_scba_result_classification(
        false,
        :max_poisson_iterations,
        :strictly_converged,
    ) == :strictly_converged
    @test QCLNEGF._check_scba_result_classification(
        false,
        :max_iterations,
        :approximate_fixed_point,
    ) == :approximate_fixed_point
    @test_throws ArgumentError QCLNEGF._check_scba_result_classification(
        true,
        :converged,
        :approximate_fixed_point,
    )
    @test_throws ArgumentError QCLNEGF._check_scba_result_classification(
        false,
        :converged,
        :strictly_converged,
    )
    @test_throws ArgumentError QCLNEGF._check_scba_result_classification(
        false,
        :max_iterations,
        :unknown,
    )

    shape = (1, 1, 1, 1)
    family = SelfEnergyFamily(
        zeros(ComplexF64, shape),
        zeros(ComplexF64, shape),
        zeros(ComplexF64, shape),
    )
    green = GreenState(
        zeros(ComplexF64, shape),
        zeros(ComplexF64, shape),
        zeros(ComplexF64, shape),
        zeros(ComplexF64, shape),
        zeros(1, 1),
        ones(1, 1),
    )
    @test_throws ArgumentError SCBAResult(
        green,
        Dict{Symbol,SelfEnergyFamily}(),
        family,
        family,
        family,
        SCBAIteration[],
        false,
        :max_iterations,
        :strictly_converged,
    )
end

end # independent suite
