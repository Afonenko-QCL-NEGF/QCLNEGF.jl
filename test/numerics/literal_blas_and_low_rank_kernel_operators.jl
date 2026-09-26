module Suite_T004
include("../support/common.jl")
include("../support/algorithm_modes.jl")

@testset "Literal, BLAS, and low-rank kernel operators" begin
    Random.seed!(0xA190)
    NE, Nk, Nb = 5, 3, 2
    K = randn(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
    G = randn(ComplexF64, NE, Nk, Nb, Nb)
    weights = rand(Nk)
    scale = 0.13
    reference = QCLNEGF._static_contraction(K, G, weights, scale)

    literal = QCLNEGF._production_kernel(
        K,
        AlgorithmOptions(contraction = :literal, kernel_build = :direct),
    )
    @test literal isa LiteralProductionKernel

    dense = QCLNEGF._production_kernel(
        K,
        AlgorithmOptions(contraction = :dense_blas, kernel_build = :direct),
    )
    @test production_static_contraction(dense, G, weights, scale; energy_chunk = 2) ≈
          reference rtol=3e-14 atol=3e-14

    full_rank = QCLNEGF._production_kernel(
        K,
        AlgorithmOptions(
            contraction = :low_rank,
            low_rank_relative_tolerance = 0.0,
            kernel_build = :direct,
        ),
    )
    @test full_rank isa LowRankProductionKernel
    @test full_rank.rank == Nk * Nb^2
    @test full_rank.relative_frobenius_residual == 0.0
    @test production_static_contraction(full_rank, G, weights, scale; energy_chunk = 2) ≈
          reference rtol=2e-13 atol=2e-13

    controlled = QCLNEGF._production_kernel(
        K,
        AlgorithmOptions(
            contraction = :low_rank,
            low_rank_relative_tolerance = 0.45,
            kernel_build = :direct,
        ),
    )
    dense_matrix = dense.matrix
    approximation = controlled.left * controlled.right_adjoint
    measured = norm(approximation - dense_matrix) / norm(dense_matrix)
    @test measured <= 0.45 + 128eps(Float64)
    @test measured ≈ controlled.relative_frobenius_residual rtol=2e-13
end

end # independent suite
