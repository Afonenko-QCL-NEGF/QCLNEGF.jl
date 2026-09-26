module Suite_T089
include("../support/common.jl")
include("../support/production_residuals.jl")

@testset "Production residual non-finite policy and contracts" begin
    ε, h, green, total, candidate = _random_residual_fixture(NE = 5, Nk = 3, Nb = 2)
    workspace = ProductionResidualWorkspace(5, 3, 2; chunk_size = 2)

    bad_candidate =
        SelfEnergyFamily(copy(candidate.Σᴿ), copy(candidate.Σˡ), copy(candidate.Σᵍ))
    bad_candidate.Σˡ[2, 1, 1, 1] = NaN + 0im
    candidate_result =
        production_residual_suite(ε, h, green, total, bad_candidate; workspace)
    @test isinf(candidate_result.r_K)
    @test isfinite(candidate_result.r_D)
    @test isfinite(candidate_result.r_A)
    @test isfinite(candidate_result.r_PSD)
    @test isfinite(candidate_result.r_caus)

    bad_total = SelfEnergyFamily(copy(total.Σᴿ), copy(total.Σˡ), copy(total.Σᵍ))
    bad_total.Σˡ[1, 2, 1, 2] = Inf + 0im
    total_result = production_residual_suite(ε, h, green, bad_total, candidate; workspace)
    @test isinf(total_result.r_A)
    @test isinf(total_result.r_PSD)
    @test isfinite(total_result.r_D)
    @test isfinite(total_result.r_caus)

    bad_retarded = SelfEnergyFamily(copy(total.Σᴿ), copy(total.Σˡ), copy(total.Σᵍ))
    bad_retarded.Σᴿ[4, 3, 2, 1] = NaN + 0im
    retarded_result =
        production_residual_suite(ε, h, green, bad_retarded, candidate; workspace)
    @test isinf(retarded_result.r_D)
    @test isinf(retarded_result.r_caus)
    @test isfinite(retarded_result.r_A)

    @test isinf(production_selfenergy_residual(total, bad_total; workspace))

    # Finite inputs may still overflow while Hermitian parts are assembled.
    # They must fail closed before entering raw LAPACK.
    overflow_A = copy(green.A)
    overflow_A[1, 1, 1, 1] = complex(floatmax(Float64), 0.0)
    overflow_green = GreenState(
        copy(green.Gᴿ),
        copy(green.Gˡ),
        copy(green.Gᵍ),
        overflow_A,
        copy(green.condition_number),
        copy(green.dyson_scale),
    )
    overflow_result =
        production_residual_suite(ε, h, overflow_green, total, candidate; workspace)
    @test isinf(overflow_result.r_PSD)

    overflow_ΣR = copy(total.Σᴿ)
    overflow_ΣR[1, 1, 1, 2] = complex(floatmax(Float64), 0.0)
    overflow_ΣR[1, 1, 2, 1] = complex(-floatmax(Float64), 0.0)
    overflow_total = SelfEnergyFamily(overflow_ΣR, copy(total.Σˡ), copy(total.Σᵍ))
    overflow_causality =
        production_residual_suite(ε, h, green, overflow_total, candidate; workspace)
    @test isinf(overflow_causality.r_caus)

    # The Dyson norm uses scaled accumulation: finite large/small operands do
    # not spuriously overflow their individual Frobenius norms and fail open.
    extreme_GR = fill(complex(1e-100, 0.0), 1, 1, 1, 1)
    extreme_green = GreenState(
        extreme_GR,
        zeros(ComplexF64, 1, 1, 1, 1),
        zeros(ComplexF64, 1, 1, 1, 1),
        zeros(ComplexF64, 1, 1, 1, 1),
        ones(1, 1),
        ones(1, 1),
    )
    extreme_family = SelfEnergyFamily(
        zeros(ComplexF64, 1, 1, 1, 1),
        zeros(ComplexF64, 1, 1, 1, 1),
        zeros(ComplexF64, 1, 1, 1, 1),
    )
    extreme_workspace = ProductionResidualWorkspace(1, 1, 1; worker_count = 1)
    extreme = production_residual_suite(
        [1e200],
        zeros(ComplexF64, 1, 1, 1),
        extreme_green,
        extreme_family,
        extreme_family;
        workspace = extreme_workspace,
    )
    @test isfinite(extreme.r_D)
    @test extreme.r_D ≈ 1.0 rtol=2e-14

    @test_throws ArgumentError ProductionResidualWorkspace(0, 3, 2)
    @test_throws ArgumentError ProductionResidualWorkspace(5, 3, 2; chunk_size = 0)
    @test_throws ArgumentError ProductionResidualWorkspace(5, 3, 2; worker_count = -1)
    wrong_workspace = ProductionResidualWorkspace(4, 3, 2)
    @test_throws DimensionMismatch production_residual_suite(
        ε,
        h,
        green,
        total,
        candidate;
        workspace = wrong_workspace,
    )
    @test_throws DimensionMismatch production_selfenergy_residual(
        total,
        candidate;
        workspace = wrong_workspace,
    )

    # Workspace validation precedes every raw LAPACK call; malformed public
    # storage must raise rather than permit an out-of-bounds Fortran write.
    short_lapack = ProductionResidualWorkspace(5, 3, 2; worker_count = 1)
    resize!(short_lapack.lapack_work[1], 1)
    @test_throws DimensionMismatch production_residual_suite(
        ε,
        h,
        green,
        total,
        candidate;
        workspace = short_lapack,
    )
    short_blocks = ProductionResidualWorkspace(5, 3, 2; worker_count = 1)
    resize!(short_blocks.spectral_numerator, 1)
    @test_throws DimensionMismatch production_selfenergy_residual(
        total,
        candidate;
        workspace = short_blocks,
    )
end

end # independent suite
