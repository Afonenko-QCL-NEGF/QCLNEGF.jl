module CavityCurrentAndSpectralOracles
include("../support/common.jl")
const BN = QCLNEGF

@testset "Cavity current equals independent coherent Landauer transmission" begin
    levels = [0.1, -0.05, 0.2]
    t1, t2, gammaL, gammaR = 0.25, 0.32, 0.07, 0.11
    onsite = [fill(ComplexF64(v), 1, 1) for v in levels]
    hopping = [fill(ComplexF64(v), 1, 1) for v in (t1, t2)]
    for energy in (-0.5, -0.1, 0.0, 0.2, 0.6), (fL, fR) in ((0.8, 0.2), (0.4, 0.4))
        sr = [fill(-0.5im*gammaL, 1, 1), zeros(ComplexF64, 1, 1), fill(-0.5im*gammaR, 1, 1)]
        sl = [fill(im*gammaL*fL, 1, 1), zeros(ComplexF64, 1, 1), fill(im*gammaR*fR, 1, 1)]
        result = BN.finite_chain_green(
            energy,
            onsite,
            hopping;
            sigma_retarded = sr,
            sigma_lesser = sl,
            centre = 2,
        )
        # Analytic inverse element G13 = t1*t2/det(D) for this three-site chain.
        a = energy-levels[1]+0.5im*gammaL
        b = energy-levels[2]
        c = energy-levels[3]+0.5im*gammaR
        determinant = a*b*c-t1^2*c-t2^2*a
        transmission = gammaL*gammaR*abs2(t1*t2/determinant)
        incoming(embedding) = real(
            only(embedding.lesser)*only(result.greater)-only(embedding.greater)*only(
                result.lesser,
            ),
        )
        @test incoming(result.minus) ≈ transmission*(fL-fR) atol=2e-14 rtol=2e-13
        @test incoming(result.plus) ≈ -transmission*(fL-fR) atol=2e-14 rtol=2e-13
        @test incoming(result.minus)+incoming(result.plus) ≈ 0 atol=3e-14
        @test 0 <= transmission <= 1+2e-14
    end
end

@testset "Finite-chain spectral integral equals analytic eigenmode weights" begin
    sites, centre, hopping, width = 5, 3, 0.3, 0.08
    levels = zeros(sites)
    onsite = [zeros(ComplexF64, 1, 1) for _ = 1:sites]
    couplings = [fill(ComplexF64(hopping), 1, 1) for _ = 1:(sites-1)]
    sr = [fill(-im*width, 1, 1) for _ = 1:sites]
    sl = [fill(0.6im*width, 1, 1) for _ = 1:sites]
    # Open uniform-chain standing waves are an independent continuum-E oracle.
    eigenvalues = [2hopping*cos(j*π/(sites+1)) for j = 1:sites]
    weights = [2/(sites+1)*sin(centre*j*π/(sites+1))^2 for j = 1:sites]
    lower, upper = -2.0, 2.0
    exact = sum(
        weights[j]*(atan((upper-eigenvalues[j])/width)-atan((lower-eigenvalues[j])/width))/π for j = 1:sites
    )
    errors = Float64[]
    for nodes in (129, 257, 513)
        energies = range(lower, upper; length = nodes)
        step = (upper-lower)/(nodes-1)
        integrated = 0.0
        for (j, energy) in enumerate(energies)
            result = BN.finite_chain_green(
                energy,
                onsite,
                couplings;
                sigma_retarded = sr,
                sigma_lesser = sl,
                centre,
            )
            spectral = -2imag(only(result.retarded))
            integrated += (j==1 || j==nodes ? 0.5 : 1.0)*step*spectral/(2π)
        end
        push!(errors, abs(integrated-exact))
    end
    @test errors[3]<errors[2]<errors[1]
    @test errors[3] < 1e-6
    @test exact < 1 # finite spectral window has an analytic tail; don't force unit area
end
end
