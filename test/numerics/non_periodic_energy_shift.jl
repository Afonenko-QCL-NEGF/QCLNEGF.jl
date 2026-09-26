module Suite_T069
include("../support/common.jl")

@testset "Non-periodic energy shift" begin
    ε = collect(0.0:0.25:1.0)
    Wplus = build_shift_matrix(ε, 0.25)
    Wminus = build_shift_matrix(ε, -0.25)
    x = reshape(ComplexF64.(ε), length(ε), 1, 1, 1)
    yplus = apply_energy_shift(Wplus, x)[:, 1, 1, 1]
    yminus = apply_energy_shift(Wminus, x)[:, 1, 1, 1]
    @test yplus[1:4] ≈ ε[2:5]
    @test yplus[end] == 0
    @test yminus[2:5] ≈ ε[1:4]
    @test yminus[1] == 0
end

end # independent suite
