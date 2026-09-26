module PSDBackwardErrorContracts
include("../support/common.jl")

@testset "PSD mixed absolute-relative backward error, independent analytic eigenvalues" begin
    @test PSD_METRIC_VERSION == 3
    @test QCLNEGF._negative_part(zeros(ComplexF64, 2, 2)) == 0
    @test QCLNEGF._negative_part(Diagonal(ComplexF64[1e-16, -1e-17])) == 0
    @test QCLNEGF._negative_part(Diagonal(ComplexF64[1.0, -1e-5])) > 1e-6
    @test QCLNEGF._negative_part(
        Diagonal(ComplexF64[1e-16, -1e-17]);
        construction_scale = 1e-16,
    ) > 0.09
    # Independent BigFloat 2×2 analytic determinant/eigenvalue witness.
    a, b, d = BigFloat(1), BigFloat("0.2"), BigFloat("0.01")
    exact = (a+d-sqrt((a-d)^2+4b^2))/2
    matrix = ComplexF64[a b; b d]
    witness = QCLNEGF._positivity_witness(:spectral, 1, 1, eigvals(Hermitian(matrix)))
    @test witness.minimum_eigenvalue ≈ Float64(exact) atol=2e-16
    @test witness.absolute_defect > witness.backward_error
    @test witness.ratio ≈
          (witness.absolute_defect-witness.backward_error)/witness.block_norm
    for scale in (1e-9, 1.0, 1e9)
        @test QCLNEGF._negative_part(scale .* matrix; construction_scale = scale) ≈
              QCLNEGF._negative_part(matrix; construction_scale = 1.0) rtol=1e-12
    end
    @test_throws ArgumentError psd_error_budget(1.0, 0)
    @test_throws ArgumentError psd_error_budget(1.0, 2; construction_scale = -1)
end
end
