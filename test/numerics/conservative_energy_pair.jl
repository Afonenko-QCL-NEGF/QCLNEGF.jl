module Suite_T108
include("../support/common.jl")

@testset "conservative energy pair" begin
    energy = [-2.0, -1.1, -0.2, 0.3, 1.8, 3.0]
    edges = vcat(first(energy), (energy[1:(end-1)]+energy[2:end])/2, last(energy))
    weights = diff(edges)
    pair = QCLNEGF.build_conservative_shift_pair(energy, weights, 0.371)
    W = Diagonal(weights)
    @test W * pair.plus ≈ pair.minus' * W atol=2e-15
    @test all(pair.plus .>= 0) && all(pair.minus .>= 0)
    @test all(sum(pair.plus; dims = 2) .<= 1+1e-14)
    @test all(sum(pair.minus; dims = 2) .<= 1+1e-14)
    @test pair.plus[end, 1] == 0 && pair.minus[1, end] == 0
    @test pair.loss_plus[end] > 0 && pair.loss_minus[1] > 0
    opposite = QCLNEGF.build_conservative_shift_pair(energy, weights, -0.371)
    @test pair.minus ≈ opposite.plus atol=2e-15
    @test QCLNEGF.build_conservative_shift_pair(energy, weights, 0).plus ≈ I
    a, b = sin.(energy), cos.(energy)
    @test dot(a, W*pair.plus*b) ≈ dot(pair.minus*a, W*b) atol=2e-15
    rescaled =
        QCLNEGF.build_conservative_shift_pair(0.3 .* energy, 0.3 .* weights, 0.3*0.371)
    @test rescaled.plus ≈ pair.plus atol=3e-15
    @test rescaled.minus ≈ pair.minus atol=3e-15
    fine_grid = collect(range(-1.5, 4.5; length = 1201))
    fine_weights = fill(6/1200, 1201)
    fine_weights[[1, end]] ./= 2
    fine_pair = QCLNEGF.build_conservative_shift_pair(fine_grid, fine_weights, 0.037)
    @test fine_pair.weighted_adjoint_residual < 1e-14

    # Independent zeroth/first moments for an interior smooth density. Zero
    # exterior values have negligible support here, so the shift must retain
    # mass and change the first moment by the signed displacement.
    nodes=collect(range(-1.0, 1.0; length = 401))
    w=fill(2/400, 401)
    w[[1, end]]./=2
    density=exp.(-((nodes .+ 0.1) ./ 0.07) .^ 2)
    mass=dot(w, density)
    moment=dot(w, nodes .* density)
    delta=0.173
    conservative=QCLNEGF.build_conservative_shift_pair(nodes, w, delta)
    for (plus, minus) in (
        (build_shift_matrix(nodes, delta), build_shift_matrix(nodes, -delta)),
        (conservative.plus, conservative.minus),
    )
        @test dot(w, plus*density) ≈ mass atol=1e-13
        @test dot(w, minus*density) ≈ mass atol=1e-13
        @test dot(w, nodes .* (plus*density)) ≈ moment-delta*mass atol=1e-13
        @test dot(w, nodes .* (minus*density)) ≈ moment+delta*mass atol=1e-13
    end
end

end # independent suite
