module Suite_T047
include("../support/common.jl")

@testset "Four-index static contraction order" begin
    NE, Nk, Nb = 2, 3, 2
    K = zeros(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
    G = zeros(ComplexF64, NE, Nk, Nb, Nb)
    w = [0.2, 0.3, 0.5]
    for m = 1:Nk, mp = 1:Nk, a = 1:Nb, c = 1:Nb, d = 1:Nb, b = 1:Nb
        K[m, mp, a, c, d, b] = m + 2mp + 3a + 5c + 7d + 11b
    end
    for e = 1:NE, m = 1:Nk, c = 1:Nb, d = 1:Nb
        G[e, m, c, d] = e + 2m + 3c + 5d
    end
    Σ = QCLNEGF._static_contraction(K, G, w)
    reference = zero(Σ)
    for e = 1:NE, m = 1:Nk, a = 1:Nb, b = 1:Nb, mp = 1:Nk, c = 1:Nb, d = 1:Nb
        reference[e, m, a, b] += w[mp] * K[m, mp, a, c, d, b] * G[e, mp, c, d]
    end
    @test Σ == reference
end

end # independent suite
