module Suite_T074
include("../support/common.jl")

@testset "Static problem without scattering" begin
    scattering = ScatteringOptions(
        LO = false,
        acoustic = false,
        impurity = false,
        IFR = false,
        alloy = false,
    )
    problem = build_problem(numerical = tutorial_numerics(), scattering = scattering)
    @test isempty(problem.kernels.enabled)
    @test size(problem.W₊ᴱᵖ) == (problem.numerical.N_E, problem.numerical.N_E)
    report = validate_problem(problem)
    @test report.passed
    @test all(iszero, problem.W₊ᴱᵖ[end, :])
    @test all(iszero, problem.W₋ᴱᵖ[1, :])
end

end # independent suite
