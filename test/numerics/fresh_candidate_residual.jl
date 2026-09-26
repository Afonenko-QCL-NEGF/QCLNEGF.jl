module FreshCandidateDysonContracts
include("../support/common.jl")

@testset "Independent resonant scalar oracle requires the fresh candidate Dyson gate" begin
    # A tiny relative change of a large self-energy can materially move a narrow
    # pole. Stored Dyson and a raw relative Sigma residual alone cannot certify
    # the Green function against the new physical map.
    sigma = 0.99 - 0.001im
    delta = 1e-9
    denominator = 1-sigma
    green = reshape(ComplexF64[inv(denominator)], 1, 1, 1, 1)
    old = reshape(ComplexF64[sigma], 1, 1, 1, 1)
    candidate = reshape(ComplexF64[sigma+delta], 1, 1, 1, 1)
    h = zeros(ComplexF64, 1, 1, 1)
    stored_residual = QCLNEGF._dyson_residual([1.0], h, old, green)
    fresh_residual = QCLNEGF._dyson_residual([1.0], h, candidate, green)
    oracle = abs(delta/denominator)/(abs((denominator-delta)/denominator)+1)
    tolerances = SolverTolerances()
    @test stored_residual < tolerances.r_D
    @test abs(delta)/abs(sigma+delta) < tolerances.r_Σ
    @test fresh_residual > tolerances.r_K
    @test fresh_residual ≈ oracle rtol=1e-6
    limits = QCLNEGF._solution_validation_limits(
        tolerances,
        Dict(:r_D_candidate=>fresh_residual),
    )
    @test limits[:r_D_candidate] == tolerances.r_K
end
end
