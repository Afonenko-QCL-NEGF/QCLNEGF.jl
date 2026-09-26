module Suite_CompactProductionShiftStorage
include("../support/common.jl")

@testset "Compact production shifts preserve the dense map and validation" begin
    BN = QCLNEGF
    for energy in
        (collect(range(-0.7, 1.1; length = 29)), [-2.0, -1.1, -0.2, 0.3, 1.8, 3.0])
        edges = vcat(first(energy), (energy[1:(end-1)] .+ energy[2:end]) ./ 2, last(energy))
        weights = diff(edges)
        window = last(energy) - first(energy)
        values = complex.(energy .^ 2, energy)
        for delta in (-2window, -window, -0.173, 0.0, 0.257, window, 2window)
            plus, minus = build_shift_plan(energy, delta), build_shift_plan(energy, -delta)
            dense_plus, dense_minus =
                build_shift_matrix(energy, delta), build_shift_matrix(energy, -delta)
            @test size(plus) == size(dense_plus)
            @test Matrix(plus) == dense_plus
            @test plus * values ≈ dense_plus * values rtol=4e-15 atol=4e-15
            @test Matrix(copy(plus)) == dense_plus
            compact_validation = BN._shift_pair_validation(
                plus,
                minus,
                energy,
                weights,
                delta,
                :nodal_linear,
            )
            dense_validation = BN._shift_pair_validation(
                dense_plus,
                dense_minus,
                energy,
                weights,
                delta,
                :nodal_linear,
            )
            for key in keys(dense_validation)
                @test getproperty(compact_validation, key) ≈
                      getproperty(dense_validation, key) atol=2e-14
            end
        end
    end
    for (energy, delta) in (([0.0], 0.0), ([0.0, Inf], 0.0), ([0.0, 1.0], Inf))
        @test_throws ArgumentError build_shift_plan(energy, delta)
    end

    energy = collect(range(-0.6, 0.6; length = 30_721))
    weights = fill(1.2 / 30_720, length(energy))
    weights[[1, end]] ./= 2
    plus, minus = build_shift_plan(energy, 0.0367), build_shift_plan(energy, -0.0367)
    @test Base.summarysize(plus) < 40length(energy)
    @test plus[end, 1] == 0.0
    @test minus[1, end] == 0.0
    # A quadratic scan/allocation on this actual finest campaign grid would
    # make this tiny validation regression impractical.
    audit = BN._shift_pair_validation(plus, minus, energy, weights, 0.0367, :nodal_linear)
    @test audit.wrap == audit.nonnegative == audit.row_mass == 0.0
    @test isfinite(audit.quadrature_adjoint)

    wrapped = copy(plus)
    wrapped.valid[end] = true
    wrapped.left[end] = wrapped.right[end] = 1
    wrapped.weight_right[end] = 0.0
    @test BN._shift_pair_validation(
        wrapped,
        minus,
        energy,
        weights,
        0.0367,
        :nodal_linear,
    ).wrap == 1.0
    negative = copy(plus)
    inside = findfirst(negative.valid)
    negative.weight_right[inside] = 1.25
    @test BN._shift_pair_validation(
        negative,
        minus,
        energy,
        weights,
        0.0367,
        :nodal_linear,
    ).nonnegative == 0.25
    malformed = copy(plus)
    malformed.left[inside] = length(energy) + 1
    @test isinf(
        BN._shift_pair_validation(malformed, minus, energy, weights, 0.0367, :nodal_linear).wrap,
    )

    empty_scattering = ScatteringOptions(
        LO = false,
        acoustic = false,
        impurity = false,
        IFR = false,
        alloy = false,
    )
    numerical = tutorial_numerics()
    direct = build_problem(numerical = numerical, scattering = empty_scattering)
    production =
        build_problem_production(numerical = numerical, scattering = empty_scattering).problem
    @test direct.W₊ᴱᵖ isa Matrix{Float64}
    @test production.W₊ᴱᵖ isa EnergyShiftPlan
    @test isconcretetype(typeof(production))
    @test validate_problem(production).passed
    for field in (:W₊ᴱᵖ, :W₋ᴱᵖ, :W₊ᴸᴼ, :W₋ᴸᴼ)
        @test Matrix(getfield(production, field)) == getfield(direct, field)
    end
    retained = retarget_problem(production; V_period = 54.0u"mV")
    @test retained.W₊ᴱᵖ isa EnergyShiftPlan
    @test retained.W₊ᴸᴼ === production.W₊ᴸᴼ
    @test validate_problem(retained).passed
    dense = retarget_problem(production; energy_shift = :dense)
    compact = retarget_problem(dense; energy_shift = :sparse_plan)
    @test dense.W₊ᴱᵖ isa Matrix{Float64}
    @test compact.W₊ᴱᵖ isa EnergyShiftPlan
    @test Matrix(compact.W₊ᴸᴼ) == dense.W₊ᴸᴼ
    @test validate_problem(compact).passed
    dense_estimate = estimate_production_memory(dense)
    compact_estimate = estimate_production_memory(compact)
    @test dense_estimate.breakdown[:dense_shift_matrices_in_problem] ==
          4numerical.N_E^2 * sizeof(Float64)
    @test compact_estimate.breakdown[:dense_shift_matrices_in_problem] == 0
    @test compact_estimate.breakdown[:compact_shift_operators_in_problem] >=
          4Base.summarysize(compact.W₊ᴱᵖ)
end
end
