module Suite_T042
include("../support/common.jl")

@testset "Energy shift validation follows the declared discretization" begin
    BN = QCLNEGF
    energy = [-2.0, -1.1, -0.2, 0.3, 1.8, 3.0]
    edges = vcat(first(energy), (energy[1:(end-1)] + energy[2:end]) / 2, last(energy))
    weights = diff(edges)
    for scale in (0.3, 1.0, 7.0), displacement in (-0.371, 0.0, 0.371, 9.0)
        ε, w, δ = scale .* energy, scale .* weights, scale * displacement
        pair = BN.build_conservative_shift_pair(ε, w, δ)
        finite_volume = BN._shift_pair_validation(
            pair.plus,
            pair.minus,
            ε,
            w,
            δ,
            :finite_volume_piecewise_constant,
        )
        @test all(value <= 1e-14 for value in values(finite_volume))
        nodal = BN._shift_pair_validation(
            build_shift_matrix(ε, δ),
            build_shift_matrix(ε, -δ),
            ε,
            w,
            δ,
            :nodal_linear,
        )
        @test all(
            getproperty(nodal, name) <= 1e-14 for
            name in (:wrap, :nonnegative, :row_mass, :cell_weights)
        )
        forward = Diagonal(w) * build_shift_matrix(ε, δ)
        backward = transpose(build_shift_matrix(ε, -δ)) * Diagonal(w)
        scale_norm = sqrt(norm(forward)^2 + norm(backward)^2)
        oracle = iszero(scale_norm) ? 0.0 : norm(forward-backward)/scale_norm
        @test nodal.quadrature_adjoint ≈ oracle atol=1e-14
    end

    # Exact failed cavity grid from the diagnostic run: an out-of-window
    # centre still receives the nonzero intersection of its boundary cell.
    ε = collect(range(-1.0, 8.0; length = 1025))
    w = fill(9 / 1024, 1025)
    w[[1, end]] ./= 2
    δ = 0.056
    pair = BN.build_conservative_shift_pair(ε, w, δ)
    @test pair.minus[7, 1] ≈ 0.12844444444443878 atol=2e-14
    @test ε[7] - δ < first(ε)
    valid = BN._shift_pair_validation(
        pair.plus,
        pair.minus,
        ε,
        w,
        δ,
        :finite_volume_piecewise_constant,
    )
    @test all(value <= 1e-14 for value in values(valid))
    @test BN._shift_pair_validation(pair.plus, pair.minus, ε, w, δ, :nodal_linear).wrap >
          0.1
    @test maximum(
        maximum(abs, view(pair.minus, i, :)) for i in eachindex(ε) if ε[i] - δ < first(ε)
    ) > 0.1

    # An actual wrap is invalid even if its opposite is set to preserve the
    # quadrature-adjoint identity. The boundary check has not been relaxed.
    wrapped_plus, wrapped_minus = copy(pair.plus), copy(pair.minus)
    wrapped_plus[end, 1] = 0.2
    wrapped_minus[1, end] = w[end] / w[1] * 0.2
    wrapped = BN._shift_pair_validation(
        wrapped_plus,
        wrapped_minus,
        ε,
        w,
        δ,
        :finite_volume_piecewise_constant,
    )
    @test wrapped.wrap == 0.2
    @test wrapped.quadrature_adjoint < 1e-14

    # A supported perturbation must still respect positivity and the weighted
    # opposite map. These are separate properties from absence of wrap.
    index = findfirst(>(0.0), pair.minus)
    broken_minus = copy(pair.minus)
    broken_minus[index] += 0.2
    broken = BN._shift_pair_validation(
        pair.plus,
        broken_minus,
        ε,
        w,
        δ,
        :finite_volume_piecewise_constant,
    )
    @test broken.wrap == 0.0
    @test broken.quadrature_adjoint > 1e-14
    broken_minus[index] = -0.2
    @test BN._shift_pair_validation(
        pair.plus,
        broken_minus,
        ε,
        w,
        δ,
        :finite_volume_piecewise_constant,
    ).nonnegative == 0.2

    # Reproduce the 481-node research grid: a finite-window nodal pair is
    # supported and nonnegative but has a measurable boundary defect.
    ε = collect(range(-0.15, 0.45; length = 481))
    w = fill(0.6/480, 481)
    w[[1, end]] ./= 2
    nodal = BN._shift_pair_validation(
        build_shift_matrix(ε, 0.0367),
        build_shift_matrix(ε, -0.0367),
        ε,
        w,
        0.0367,
        :nodal_linear,
    )
    @test nodal.quadrature_adjoint ≈ 0.0235636814813 rtol=1e-9
    @test nodal.wrap <= 1e-14
end

end # independent suite
