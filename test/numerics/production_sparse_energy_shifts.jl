module Suite_T076
include("../support/common.jl")
include("../support/production_backend.jl")

@testset "Production sparse energy shifts" begin
    ε = collect(range(-0.7, 1.1; length = 29))
    X = reshape(ComplexF64.(1:(29*3*2*2)), 29, 3, 2, 2)
    for δ in (-0.173, 0.0, 0.257, ε[end] - ε[1])
        dense = build_shift_matrix(ε, δ)
        plan = build_shift_plan(ε, δ)
        @test isapprox(
            apply_energy_shift(plan, X),
            apply_energy_shift(dense, X);
            rtol = 4e-15,
            atol = 4e-15,
        )
        destination = similar(X)
        @test apply_energy_shift!(destination, plan, X) === destination
        @test isapprox(
            destination,
            apply_energy_shift(dense, X);
            rtol = 4e-15,
            atol = 4e-15,
        )
    end
end

end # independent suite
