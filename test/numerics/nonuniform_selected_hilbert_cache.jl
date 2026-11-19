module NonuniformSelectedHilbertCache
include("../support/common.jl")

# Internal operator fixture, not a public nonuniform solver/grid contract.
function cache_fixture()
    original = tutorial_numerics()
    changes = (; N_b = 2, N_k = 2, N_E = 5, E_min = -0.2u"eV", E_max = 0.3u"eV")
    numerical = NumericalParameters(;
        (name => hasproperty(changes, name) ? getproperty(changes, name) :
                 getfield(original, name) for name in fieldnames(NumericalParameters))...,
    )
    uniform = build_problem(;
        numerical,
        scattering = ScatteringOptions(
            LO = false, acoustic = false, impurity = false, IFR = false, alloy = false,
        ),
    )
    energy = [-2.0, -0.7, 0.0, 0.4, 3.0]
    weights = vcat(diff(energy)[1]/2, (energy[3:end] - energy[1:end-2])/2,
                   diff(energy)[end]/2)
    g = uniform.grids
    grids = QCLNEGF.QCLDomain.ModelGrids(g.x, g.wˣ, energy, weights, trues(5), g.κ, g.wᵏ,
                       g.qᶻ, g.wᑫᶻ, g.φ)
    sp = QCLNEGF._scaled_physics(uniform.physical, uniform.scales)
    shifts = map(δ -> build_shift_plan(energy, δ),
                 (+sp.Eᵖ, -sp.Eᵖ, +sp.ħωᴸᴼ, -sp.ħωᴸᴼ))
    nonuniform = NEGFProblem(
        uniform.physical, numerical, uniform.scattering, uniform.scales, grids,
        uniform.profiles, uniform.basis, uniform.kernels, shifts...,
        :nodal_linear, uniform.models,
    )
    return uniform, nonuniform
end

function broadening_tensor(values)
    block = ComplexF64[1 0.2+0.1im; 0.2-0.1im 2]
    result = zeros(ComplexF64, length(values), 2, 2, 2)
    for e in eachindex(values), k in 1:2
        result[e, k, :, :] .= k * values[e] .* block
    end
    return result
end

# Closed form of the triangular hat's finite-support PV integral, divided by
# 2π. Unlike the production matrix, this uses the three slope jumps directly.
xlogabs(x) = iszero(x) ? 0.0 : x * log(abs(x))
triangle_pv(E, a, c, b) = (
    xlogabs(E-a)/(c-a) - (1/(c-a)+1/(b-c))*xlogabs(E-c) + xlogabs(E-b)/(b-c)
)/(2π)

@testset "Nonuniform selected Hilbert cache and retarded operator" begin
    uniform, problem = cache_fixture()
    energy = problem.grids.ε
    fft_estimate = estimate_production_memory(problem)
    direct_estimate = estimate_production_memory(problem;
        options = ProductionOptions(algorithms = AlgorithmOptions(hilbert = :direct)))
    product_estimate = estimate_production_memory(problem;
        options = ProductionOptions(algorithms = AlgorithmOptions(hilbert = :product_integration)))
    # Only the selected persistent Hilbert storage contributes to the RAM guard.
    @test direct_estimate.breakdown[:fft_workspace] == 0
    @test product_estimate.breakdown[:fft_workspace] == 0
    @test fft_estimate.breakdown[:fft_workspace] > 0
    @test fft_estimate.peak_bytes - direct_estimate.peak_bytes ==
          fft_estimate.breakdown[:fft_workspace]
    dense_bytes = 5^2*sizeof(Float64)
    @test product_estimate.resident_bytes - direct_estimate.resident_bytes == dense_bytes
    @test product_estimate.peak_bytes - direct_estimate.peak_bytes == 2dense_bytes
    # At large NE the two simultaneously live H arrays dominate. This is
    # arithmetic planning only; do not allocate a large dense operator here.
    n = problem.numerical
    large = NumericalParameters(;
        (name => name === :N_E ? 4097 : getfield(n, name)
         for name in fieldnames(NumericalParameters))...,
    )
    large_estimate = estimate_production_memory(large, 0;
        dense_mechanism_count = 0,
        options = ProductionOptions(algorithms = AlgorithmOptions(hilbert = :product_integration)),
    )
    @test large_estimate.peak_bytes >= 2*4097^2*sizeof(Float64)
    @test large_estimate.breakdown[:product_hilbert_verification_temporary] ==
          large_estimate.breakdown[:product_hilbert_operator] == 4097^2*sizeof(Float64)

    for mode in (:direct, :product_integration)
        estimate = mode === :direct ? direct_estimate : product_estimate
        options = ProductionOptions(;
            algorithms = AlgorithmOptions(hilbert = mode),
            memory_budget_bytes = estimate.peak_bytes,
        )
        cache = build_production_cache(problem; options)
        cache = QCLNEGF._validated_production_cache(cache, problem, options)
        too_small = QCLNEGF.with_production_options(options;
            memory_budget_bytes = estimate.peak_bytes-1)
        @test_throws ArgumentError build_production_cache(problem; options = too_small)
        @test_throws ArgumentError QCLNEGF._validated_production_cache(cache, problem, too_small)
        if mode === :product_integration
            # The previous admission threshold omitted the simultaneously live
            # reversed H. Both construction and cache reuse must reject it.
            missing_reverse_H = QCLNEGF.with_production_options(options;
                memory_budget_bytes = estimate.peak_bytes-dense_bytes)
            @test_throws ArgumentError build_production_cache(problem; options = missing_reverse_H)
            @test_throws ArgumentError QCLNEGF._validated_production_cache(
                cache, problem, missing_reverse_H)
        end
        @test_throws ArgumentError QCLNEGF._validated_production_cache(cache, problem,
            QCLNEGF.with_production_options(options; algorithms = AlgorithmOptions(hilbert = :fft)))
        if mode === :direct
            # At E=0, Γ(x)=x² cancels the singularity: the integrand is -x.
            # Nonuniform trapezoids integrate it exactly on [-2,3], giving
            # PV ∫ x²/(2π*(0-x)) dx = (a²-b²)/(4π). This checks the retained
            # omitted-node rule at a removable singularity, not general
            # continuum consistency of direct PV on arbitrary grids.
            values = energy .^ 2
            exact = (first(energy)^2-last(energy)^2)/(4π)
        else
            # The sampled triangle is exactly piecewise linear, zero at the
            # endpoints, so the production exterior taper adds no contribution.
            values = [E <= 0 ? (E+2)/2 : (3-E)/3 for E in energy]
            exact = [triangle_pv(E, -2.0, 0.0, 3.0) for E in energy]
        end
        Γ = broadening_tensor(values)
        retarded, actual_Γ, roundoff = QCLNEGF._retarded_self_energy_selected(
            0.3im .* Γ, -0.7im .* Γ, problem, cache, options; return_roundoff = true,
        )
        @test actual_Γ ≈ Γ atol=2e-15 rtol=2e-15
        @test isfinite(roundoff) && roundoff < 2e-14
        @test all(isfinite, retarded)
        if mode === :direct
            expected = broadening_tensor(fill(exact, 5)) - 0.5im .* Γ
            @test retarded[3, :, :, :] ≈ expected[3, :, :, :] atol=3e-14 rtol=3e-14
        else
            @test retarded ≈ broadening_tensor(exact) - 0.5im .* Γ atol=3e-14 rtol=3e-14
        end
    end

    # Default FFT selection retains its uniform-grid and cache-option guards.
    @test_throws ArgumentError build_production_cache(problem)
    options = ProductionOptions()
    cache = build_production_cache(uniform; options)
    QCLNEGF._validated_production_cache(cache, uniform, options)
    @test_throws ArgumentError QCLNEGF._validated_production_cache(
        cache, uniform, QCLNEGF.with_production_options(options; hilbert_columns = 1),
    )
end
end
