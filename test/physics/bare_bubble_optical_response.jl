module Suite_T083
include("../support/common.jl")

@testset "Bare-bubble optical response" begin
    εsmall = collect(range(-1.0, 1.0; length = 9))
    left, θ, valid = QCLNEGF._optical_shift_plan(εsmall, 0.375)
    @test valid == BitVector([true, true, true, true, true, true, true, false, false])
    for e in findall(valid)
        reconstructed = (1 - θ[e]) * εsmall[left[e]] + θ[e] * εsmall[left[e]+1]
        @test reconstructed ≈ εsmall[e] + 0.375 atol=5e-16
    end
    @test_throws ArgumentError QCLNEGF._optical_shift_plan(εsmall, -0.1)

    # Exactly solvable two-level equilibrium fixture.  It checks the response
    # sign without inserting a phenomenological gain formula: lower-state
    # occupation must absorb, inversion must amplify, and equal occupations
    # must have no dissipative response.
    ε = collect(range(-2.0, 3.0; length = 2001))
    Δε = ε[2] - ε[1]
    wᴱ = fill(Δε, length(ε))
    wᴱ[[1, end]] ./= 2
    wᵏ = [1.0]
    Z = ComplexF64[0 1; 1 0]
    E₁, E₂, γ = 0.2, 0.8, 0.02

    function manufactured_green(f₁, f₂)
        Gᴿ = zeros(ComplexF64, length(ε), 1, 2, 2)
        Gˡ = similar(Gᴿ)
        fill!(Gˡ, 0)
        for e in eachindex(ε)
            Gᴿ[e, 1, 1, 1] = inv(ε[e] - E₁ + im * γ)
            Gᴿ[e, 1, 2, 2] = inv(ε[e] - E₂ + im * γ)
            for a = 1:2
                Aaa = -2imag(Gᴿ[e, 1, a, a])
                Gˡ[e, 1, a, a] = im * (a == 1 ? f₁ : f₂) * Aaa
            end
        end
        return Gᴿ, Gˡ
    end

    Gᴿabs, Gˡabs = manufactured_green(1.0, 0.0)
    Sabs, _ = QCLNEGF._bare_bubble_sum(ε, wᴱ, wᵏ, Z, Gᴿabs, Gˡabs, E₂ - E₁)
    ZGᴿZ = QCLNEGF._dipole_sandwich_blocks(Z, Gᴿabs)
    ZGˡZ = QCLNEGF._dipole_sandwich_blocks(Z, Gˡabs)
    Sabs_cached, _ =
        QCLNEGF._bare_bubble_sum_cached(ε, wᴱ, wᵏ, ZGᴿZ, ZGˡZ, Gᴿabs, Gˡabs, E₂ - E₁)
    @test Sabs_cached ≈ Sabs rtol=2e-14
    Gᴿgain, Gˡgain = manufactured_green(0.0, 1.0)
    Sgain, _ = QCLNEGF._bare_bubble_sum(ε, wᴱ, wᵏ, Z, Gᴿgain, Gˡgain, E₂ - E₁)
    Gᴿflat, Gˡflat = manufactured_green(0.5, 0.5)
    Sflat, _ = QCLNEGF._bare_bubble_sum(ε, wᴱ, wᵏ, Z, Gᴿflat, Gˡflat, E₂ - E₁)

    # All omitted physical prefactors are positive except the product of the
    # electron charge and -i from the density, giving χ ∝ iS.
    @test imag(im * Sabs) > 0       # passive medium: Im χ > 0
    @test imag(im * Sgain) < 0      # inversion: Im χ < 0, hence g > 0
    @test abs(imag(im * Sflat)) ≤ 1e-10 * abs(imag(im * Sabs))
    @test Sgain ≈ -Sabs rtol=1e-5
end

end # independent suite
