module Suite_T078
include("../support/common.jl")
include("../support/production_backend.jl")

@testset "FFT linear Hilbert map equals direct principal value" begin
    Random.seed!(0xFF71)
    NE, Nk, Nb = 37, 2, 2
    ε = collect(range(-1.3, 0.9; length = NE))
    Δε = ε[2] - ε[1]
    wᴱ = fill(Δε, NE)
    wᴱ[[1, end]] ./= 2
    Γ = randn(ComplexF64, NE, Nk, Nb, Nb)
    direct = direct_hilbert_transform(Γ, ε, wᴱ)
    transformed, rround =
        fft_hilbert_transform(Γ, ε, wᴱ; column_chunk = 3, return_roundoff = true)
    @test transformed ≈ direct rtol=2e-13 atol=2e-13
    @test isfinite(rround)
    @test rround < 1e-12

    Σˡ = randn(ComplexF64, NE, Nk, Nb, Nb)
    Σᵍ = randn(ComplexF64, NE, Nk, Nb, Nb)
    ΣRdirect, Γdirect = retarded_self_energy(
        Σˡ,
        Σᵍ,
        ModelGrids(
            Float64[],
            Float64[],
            ε,
            wᴱ,
            trues(NE),
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
        ),
    )
    ΣRfft, Γfft = retarded_self_energy_fft(
        Σˡ,
        Σᵍ,
        ModelGrids(
            Float64[],
            Float64[],
            ε,
            wᴱ,
            trues(NE),
            Float64[],
            Float64[],
            Float64[],
            Float64[],
            Float64[],
        );
        column_chunk = 5,
    )
    @test Γfft == Γdirect
    @test ΣRfft ≈ ΣRdirect rtol=2e-13 atol=2e-13
end

end # independent suite
