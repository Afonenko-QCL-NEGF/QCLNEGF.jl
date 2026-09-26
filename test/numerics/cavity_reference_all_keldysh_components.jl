module Suite_T109
include("../support/common.jl")

@testset "cavity reference all Keldysh components" begin
    count, Nb = 5, 2
    onsite = [ComplexF64[0.13j 0.1im; -0.1im 0.4-0.07j] for j = 1:count]
    couplings = [ComplexF64[-0.31 0.04im; 0.06 -0.2] for _ = 1:(count-1)]
    Γ = [Matrix(Diagonal([0.07+0.01j, 0.13])) for j = 1:count]
    sigmaR = [-0.5im .* g for g in Γ]
    sigmaL = [im .* (0.1+0.12j) .* Γ[j] for j = 1:count]
    result = QCLNEGF.finite_chain_green(
        0.21,
        onsite,
        couplings;
        sigma_retarded = sigmaR,
        sigma_lesser = sigmaL,
        centre = 3,
        return_full = true,
    )
    for component in (:retarded, :lesser, :greater)
        @test getproperty(result, component) ≈ getproperty(result.full, component)[5:6, 5:6] rtol=3e-13 atol=3e-13
    end
    @test result.greater - result.lesser ≈ result.retarded - result.retarded' atol=2e-13
    @test minimum(eigvals(Hermitian(im*(result.retarded-result.retarded')))) >= 0
    # R/< /> Green functions all carry inverse-energy units, while every
    # embedding carries energy units; no arbitrary scalar scale can change
    # dimensioned observables.
    scale = 0.37
    rescaled = QCLNEGF.finite_chain_green(
        scale*0.21,
        [scale .* h for h in onsite],
        [scale .* t for t in couplings];
        sigma_retarded = [scale .* v for v in sigmaR],
        sigma_lesser = [scale .* v for v in sigmaL],
        centre = 3,
    )
    for component in (:retarded, :lesser, :greater)
        @test scale .* getproperty(rescaled, component) ≈ getproperty(result, component) rtol=5e-13
        @test getproperty(rescaled.plus, component) ./ scale ≈
              getproperty(result.plus, component) rtol=5e-13
    end

    hopping, broadening, f = 0.7, 0.12, 0.37
    for energy in (-2.0, 0.0, 0.5, 2.0)
        z = energy + im*broadening
        surface = QCLNEGF.chain_surface_green(z, hopping)
        @test surface ≈ inv(z-hopping^2*surface) rtol=2e-14
        @test imag(surface) < 0
        N = 201
        chain = QCLNEGF.finite_chain_green(
            energy,
            [zeros(ComplexF64, 1, 1) for _ = 1:N],
            [fill(ComplexF64(hopping), 1, 1) for _ = 1:(N-1)];
            sigma_retarded = [fill(-im*broadening, 1, 1) for _ = 1:N],
            sigma_lesser = [fill(2im*broadening*f, 1, 1) for _ = 1:N],
        )
        exact = inv(z-2hopping^2*surface)
        @test only(chain.retarded) ≈ exact rtol=2e-6
        @test only(chain.lesser) ≈ -2im*f*imag(exact) rtol=2e-6
        @test only(chain.greater) ≈ 2im*(1-f)*imag(exact) rtol=2e-6
    end
end

end # independent suite
