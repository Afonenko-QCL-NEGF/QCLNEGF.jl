using LinearAlgebra
using QCLNEGF
using SparseArrays
using Test

@testset "Projected spectrum uses an explicit solve and observation domain" begin
    # H is the symmetric two-site operator [2 1; 1 4] eV in coordinates
    # whitened by sqrt(M), M = diag(1, 4). Its roots are 3 ± sqrt(2) eV.
    H = [2.0 2.0; 0.5 4.0]
    W = [1.0 1.0; 0.0 1.0]
    metric = [1.0, 4.0]
    full = projected_spectrum(H, W, metric)
    cropped = projected_spectrum(H, W, metric; observation_rows = [2])
    @test full.energies_eV ≈ [3 - sqrt(2), 3 + sqrt(2)] atol = 1e-12 rtol = 0
    @test full.amplitudes' * Diagonal(metric) * full.amplitudes ≈ Matrix{Float64}(I, 2, 2)
    @test projected_spectrum(sparse(H), sparse(W), metric).energies_eV ≈ full.energies_eV
    @test cropped.energies_eV == full.energies_eV
    @test abs2.(cropped.amplitudes) ≈ abs2.(full.amplitudes[2:2, :])
    @test sum(cropped.metric .* abs2.(cropped.amplitudes[:, 1])) < 1

    # A one-column Ritz problem has the Rayleigh quotient 22/5 eV.
    one_mode = reshape([1.0, 1.0], 2, 1)
    displayed = projected_spectrum(H, one_mode, metric; observation_rows = [1])
    restricted = projected_spectrum(H, one_mode, metric; solve_rows = [1])
    @test only(displayed.energies_eV) ≈ 22 / 5 atol = 1e-12 rtol = 0
    @test only(restricted.energies_eV) ≈ 2.0 atol = 1e-12 rtol = 0
    @test abs2(only(displayed.amplitudes)) ≈ 0.2
    @test abs2(only(restricted.amplitudes)) ≈ 1.0
    @test_throws ArgumentError projected_spectrum(H, W, ones(2))
    @test_throws ArgumentError projected_spectrum(H, [1.0 1.0; 0.0 0.0], metric)
    @test_throws ArgumentError projected_spectrum(H, one_mode, metric;
        solve_rows = [1], observation_rows = [2])
    @test_throws ArgumentError projected_spectrum(H, one_mode, metric; observation_rows = [true])
    @test_throws DimensionMismatch projected_spectrum(H, W, ones(3))
end
