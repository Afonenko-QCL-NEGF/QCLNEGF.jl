module Suite_T110
include("../support/common.jl")

@testset "forward angular quadrature" begin
    # Known continuous angular average: <1/(a-b cos φ)> = 1/sqrt(a²-b²).
    # With q²=2-2cosφ this is sharply forward-peaked but nonsingular.
    a, b = 2.0004, 2.0
    evaluator(q) = fill(ComplexF64(inv(a-b+q^2)), 1, 1, 1, 1)
    quadrature =
        QCLNEGF.adaptive_angular_average(evaluator, 1.0, 1.0; relative_tolerance = 1e-7)
    @test quadrature.accepted
    @test real(only(quadrature.value)) ≈ inv(sqrt(a*a-b*b)) rtol=2e-7
    @test quadrature.estimated_relative_error > 0
    @test quadrature.evaluations > 16
    failed = QCLNEGF.adaptive_angular_average(
        evaluator,
        1.0,
        1.0;
        relative_tolerance = 1e-14,
        maximum_depth = 1,
        initial_panels = 1,
    )
    @test !failed.accepted
end

end # independent suite
