module Suite_T070
include("../support/common.jl")

@testset "Scalar Dyson and spectral identity" begin
    ε = collect(range(-1.0, 1.0; length = 2001))
    E★ = 0.17
    Γ = 0.08
    h = reshape(ComplexF64[E★], 1, 1, 1)
    Σᴿ = fill(-0.5im * Γ, length(ε), 1, 1, 1)
    Gᴿ, κD = retarded_green(ε, h, Σᴿ)
    A = spectral_function(Gᴿ)
    exact = @. Γ / ((ε - E★)^2 + (Γ / 2)^2)
    @test real.(A[:, 1, 1, 1]) ≈ exact rtol=2e-13
    @test all(κD .≈ 1)

    f = @. inv(exp((ε - 0.0) / 0.1) + 1)
    Gˡ = similar(Gᴿ)
    Gˡ[:, 1, 1, 1] .= im .* f .* A[:, 1, 1, 1]
    Gᵍ = greater_green(Gᴿ, Gˡ)
    @test maximum(abs, Gᵍ .- Gˡ .- (Gᴿ .- conj.(Gᴿ))) < 1e-12
end

end # independent suite
