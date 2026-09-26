module Suite_T020
include("../support/common.jl")

@testset "BenDaniel-Duke direct matrix" begin
    N = 12
    Δx = 0.2
    Cbar = 0.5
    mass = 0.1
    H = build_bdd_hamiltonian(zeros(N), fill(mass, N), Δx, Cbar)
    @test H ≈ H'
    @test all(iszero(H[i, j]) for i = 1:N, j = 1:N if abs(i - j) > 1)
    t = Cbar / (mass * Δx^2)
    exact = [2t * (1 - cos(j * π / (N + 1))) for j = 1:N]
    @test eigvals(Hermitian(H)) ≈ exact rtol=1e-12
end

end # independent suite
