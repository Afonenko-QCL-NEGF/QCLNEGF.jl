module Suite_T005
include("../support/common.jl")
include("../support/algorithm_modes.jl")

@testset "Physical approximation operators are explicit" begin
    Random.seed!(0x50485953)
    NE, Nk, Nb = 4, 3, 2
    family = SelfEnergyFamily(
        randn(ComplexF64, NE, Nk, Nb, Nb),
        randn(ComplexF64, NE, Nk, Nb, Nb),
        randn(ComplexF64, NE, Nk, Nb, Nb),
    )
    diagonal = QCLNEGF._copy_family(family)
    QCLNEGF._diagonalize_blocks!(diagonal.Σᴿ)
    @test all(iszero, diagonal.Σᴿ[:, :, 1, 2])
    @test all(iszero, diagonal.Σᴿ[:, :, 2, 1])
    @test diagonal.Σᴿ[:, :, 1, 1] == family.Σᴿ[:, :, 1, 1]

    weights = [0.2, 0.3, 0.5]
    averaged = copy(family.Σˡ)
    reference =
        dropdims(sum(averaged .* reshape(weights, 1, Nk, 1, 1); dims = 2); dims = 2) ./
        sum(weights)
    QCLNEGF._average_transverse_momentum!(averaged, weights)
    for m = 1:Nk
        @test averaged[:, m, :, :] ≈ reference rtol=8eps(Float64)
    end
end

end # independent suite
