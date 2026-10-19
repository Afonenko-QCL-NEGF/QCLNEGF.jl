module KeldyshLiteralMatrixOracle
include("../support/common.jl")

@testset "Literal noncommuting Keldysh product" begin
    # O08-01: stationary homogeneous functional, not a Dyson/SCBA solution.
    # Gaussian-integer arithmetic is exact in binary64; no quadrature weights.
    GR = ComplexF64[1-im 1; 0 2-im]
    C = ComplexF64[2 im; -im 1]
    lesser = im .* C
    expected = ComplexF64[7im -4+3im; 4+3im 5im]
    @test GR*C != C*GR
    actual = keldysh_green(reshape(GR, 1, 1, 2, 2), reshape(lesser, 1, 1, 2, 2))
    block = actual[1, 1, :, :]
    @test block == expected
    @test block + block' == zeros(ComplexF64, 2, 2)
    @test norm(block-expected) == 0
    B = -im .* block
    @test B == ComplexF64[7 3+4im; 3-4im 5]
    @test ishermitian(C) && ishermitian(B)
    @test C[1, 1] == 2 && C[1, 1]*C[2, 2]-C[1, 2]*C[2, 1] == 1 && tr(C) == 3
    @test B[1, 1] == 7 && B[1, 1]*B[2, 2]-B[1, 2]*B[2, 1] == 10 && tr(B) == 12
    # Positive principal minors certify positive definiteness for these 2x2
    # Hermitian matrices. PSD alone would miss a positive multiplier error.
    for mutated in (-GR*lesser*GR', 2GR*lesser*GR', GR*lesser*GR, lesser*GR*GR')
        @test mutated != expected
        @test norm(mutated-expected) > 0
    end
end
end
