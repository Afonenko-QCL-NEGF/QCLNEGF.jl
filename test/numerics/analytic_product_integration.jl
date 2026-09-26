module Suite_T107
include("../support/common.jl")

@testset "analytic product integration" begin
    e = [0.0, 1.0, 2.0]
    Γ = reshape(ComplexF64[0, 1, 0], 3, 1, 1, 1)
    result, roundoff = QCLNEGF.product_integration_hilbert_transform(
        Γ,
        e,
        [0.5, 1.0, 0.5];
        return_roundoff = true,
    )
    @test vec(real.(result)) ≈ [-log(2)/π, 0, log(2)/π] atol=2e-15
    @test roundoff < 1e-14
    nonuniform = [-2.0, -0.1, 3.0]
    answer =
        QCLNEGF.product_integration_hilbert_transform(Γ, nonuniform, [0.95, 2.5, 1.55])
    @test real(answer[2]) ≈ log(1.9/3.1)/(2π) atol=2e-15
    errors = Float64[]
    for N in (401, 801)
        energy = collect(range(-20.0, 20.0; length = N))
        broadening = reshape(ComplexF64.(2 ./ (1 .+ energy .^ 2)), N, 1, 1, 1)
        actual = QCLNEGF.product_integration_hilbert_transform(
            broadening,
            energy,
            fill(40/(N-1), N),
        )
        interior = abs.(energy) .<= 3
        exact = energy ./ (1 .+ energy .^ 2)
        push!(errors, maximum(abs, real.(vec(actual))[interior] - exact[interior]))
    end
    @test errors[2] < errors[1]/2
    @test errors[2] < 5e-4
    @test_throws ArgumentError QCLNEGF.product_integration_hilbert_matrix([0.0, 0.0, 1.0])
    # Energy-unit and reference-energy changes do not alter the PV operator.
    @test QCLNEGF.product_integration_hilbert_matrix(3.7 .* nonuniform .+ 0.41) ≈
          QCLNEGF.product_integration_hilbert_matrix(nonuniform) atol=2e-15
end

end # independent suite
