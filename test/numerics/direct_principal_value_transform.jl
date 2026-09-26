module Suite_T049
include("../support/common.jl")

@testset "Direct principal-value transform" begin
    ε = collect(-2.0:0.05:2.0)
    w = fill(0.05, length(ε))
    w[[1, end]] ./= 2
    Γ = zeros(ComplexF64, length(ε), 1, 1, 1)
    Γ[:, 1, 1, 1] .= exp.(-ε .^ 2)
    Λ = direct_hilbert_transform(Γ, ε, w)
    centre = findfirst(iszero, ε)
    @test abs(Λ[centre, 1, 1, 1]) < 1e-14
    @test real(Λ[centre+5, 1, 1, 1]) ≈ -real(Λ[centre-5, 1, 1, 1]) rtol=1e-12
end

end # independent suite
