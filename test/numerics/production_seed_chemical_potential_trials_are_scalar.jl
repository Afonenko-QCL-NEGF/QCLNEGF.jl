module Suite_T092
include("../support/common.jl")

@testset "Production seed chemical-potential trials are scalar" begin
    scattering = ScatteringOptions(
        LO = false,
        acoustic = false,
        impurity = false,
        IFR = false,
        alloy = false,
    )
    problem = build_problem(numerical = tutorial_numerics(), scattering = scattering)
    n = problem.numerical
    h = project_hamiltonians(problem, zeros(n.N_z))
    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    zeroΣ = zeros(ComplexF64, shape)
    η = scaled_value(n.η_seed, problem.scales, :energy)
    Gᴿ, _, _ = QCLNEGF._retarded_green_production(problem.grids.ε, h, zeroΣ; η = η)
    q = QCLNEGF._seed_number_weights_production(Gᴿ, problem.grids, problem.physical.g_s)
    @test size(q) == (n.N_E,)
    @test eltype(q) === Float64
    @test all(isfinite, q)

    kBT =
        QCLNEGF._electronvolts(QCLNEGF.CODATA.kᴮₑᵥ * problem.physical.Tᴸ) /
        problem.scales.E₀_eV
    μ = (problem.grids.ε[1] + problem.grids.ε[end]) / 2
    weighted_number = QCLNEGF._seed_number_from_weights(μ, problem.grids.ε, q, kBT)

    A = spectral_function(Gᴿ)
    Gˡ = zeros(ComplexF64, shape)
    for e = 1:n.N_E
        Gˡ[e, :, :, :] .= im * QCLNEGF._fermi(problem.grids.ε[e], μ, kBT) .* A[e, :, :, :]
    end
    direct_number = number_functional(Gˡ, problem.grids, problem.physical.g_s)
    @test weighted_number ≈ direct_number rtol=16eps(Float64) atol=16eps(Float64)

    # After compilation, one trial evaluates only the N_E-vector reduction.
    # In particular it cannot allocate the former full complex lesser field.
    QCLNEGF._seed_number_from_weights(μ, problem.grids.ε, q, kBT)
    trial_allocations =
        @allocated QCLNEGF._seed_number_from_weights(μ, problem.grids.ε, q, kBT)
    full_lesser_bytes = sizeof(ComplexF64) * prod(shape)
    @test trial_allocations < min(4096, full_lesser_bytes ÷ 4)

    @test_throws DimensionMismatch QCLNEGF._seed_number_from_weights(
        μ,
        problem.grids.ε,
        q[1:(end-1)],
        kBT,
    )
    @test_throws DimensionMismatch QCLNEGF._seed_lesser_production!(
        zeros(ComplexF64, n.N_E, n.N_k, n.N_b, n.N_b - 1),
        Gᴿ,
        problem.grids.ε,
        μ,
        kBT,
    )
end

end # independent suite
