module Suite_T077
include("../support/common.jl")
include("../support/production_backend.jl")

@testset "Production contractions equal literal six-index map" begin
    Random.seed!(0xB05C0)
    NE, Nk, Nb = 7, 4, 2
    G = randn(ComplexF64, NE, Nk, Nb, Nb)
    K = randn(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
    wᵏ = rand(Nk)
    qᴷ = 0.037
    direct = QCLNEGF._static_contraction(K, G, wᵏ, qᴷ)
    optimized = production_static_contraction(K, G, wᵏ, qᴷ; energy_chunk = 3)
    @test optimized ≈ direct rtol=2e-14 atol=2e-14

    block = randn(ComplexF64, Nb, Nb, Nb, Nb)
    Kindependent = zeros(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
    for m = 1:Nk, mp = 1:Nk
        Kindependent[m, mp, :, :, :, :] .= block
    end
    operator = QCLNEGF._production_kernel(Kindependent)
    @test operator isa QCLNEGF.MomentumIndependentProductionKernel
    reduced = production_static_contraction(operator, G, wᵏ, qᴷ; energy_chunk = 2)
    reference = QCLNEGF._static_contraction(Kindependent, G, wᵏ, qᴷ)
    @test reduced ≈ reference rtol=2e-14 atol=2e-14
end

end # independent suite
