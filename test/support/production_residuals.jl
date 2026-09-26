using Random

function _random_residual_fixture(; NE = 11, Nk = 5, Nb = 3, seed = 0xB05C0123)
    Random.seed!(seed)
    ε = collect(range(-0.9, 1.2; length = NE))
    h = zeros(ComplexF64, Nk, Nb, Nb)
    for m = 1:Nk
        block = randn(ComplexF64, Nb, Nb)
        h[m, :, :] .= (block + block') / 7
    end
    Gᴿ = randn(ComplexF64, NE, Nk, Nb, Nb) ./ 3
    Gˡ = randn(ComplexF64, NE, Nk, Nb, Nb) ./ 9
    Gᵍ = randn(ComplexF64, NE, Nk, Nb, Nb) ./ 8
    A = randn(ComplexF64, NE, Nk, Nb, Nb) ./ 4
    green = GreenState(Gᴿ, Gˡ, Gᵍ, A, ones(NE, Nk), ones(NE, Nk))
    family(scale) = SelfEnergyFamily(
        scale .* randn(ComplexF64, NE, Nk, Nb, Nb),
        scale .* randn(ComplexF64, NE, Nk, Nb, Nb),
        scale .* randn(ComplexF64, NE, Nk, Nb, Nb),
    )
    return ε, h, green, family(0.13), family(0.11)
end
