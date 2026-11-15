module UniformProductIntegrationFFT
include("../support/common.jl")
using Random
const N = QCLNEGF.QCLNumerics

# Positive-series oracle from even moments of the unit triangular hat.
# At n>=2, the remainder after m=159 is bounded by the first omitted term
# divided by (1-1/n²), below 1e-96 at the worst n=2. No cancelling logs.
function coefficient_oracle(n)
    n == 0 && return 0.0
    abs(n) == 1 && return sign(n)*log(2)/π
    return setprecision(256) do
        x = BigFloat(abs(n))
        value = sum(2/(BigFloat(2m+1)*(2m+2)*x^(2m+1)) for m in 0:159)
        Float64(sign(n)*value/(2big(π)))
    end
end
F(x) = iszero(x) ? zero(x) : x*log(abs(x))
function hat_oracle(E, a, c, b)
    return setprecision(256) do
        e, l, p, r = BigFloat.((E, a, c, b))
        Float64((F(e-l)/(p-l) - (1/(p-l)+1/(r-p))*F(e-p) + F(e-r)/(r-p))/(2big(π)))
    end
end

function pi_apply(Γ, energy; chunk = 3, workers = 1, weights = ones(length(energy)))
    plan = N._ProductionPIFFTHilbertPlan(energy, prod(size(Γ)[2:end]);
        column_chunk = chunk, worker_count = workers)
    result, roundoff = N._production_pi_fft_hilbert_transform(Γ, energy, weights, plan)
    return result, roundoff, plan
end

@testset "Uniform PI FFT independent coefficients, impulses and edge ramps" begin
    for lag in (0, 1, 2, 3, 7, 32, 1024, 30720, 2^30)
        coefficient = N._product_integration_lag_coefficient(lag)
        @test coefficient ≈ coefficient_oracle(lag) rtol=16eps(Float64) atol=0
        @test N._product_integration_lag_coefficient(-lag) == -coefficient
    end
    for nodes in (2, 3, 8, 17, 32, 41)
        energy = collect(range(-1.3, 0.9; length = nodes))
        impulses = reshape(Matrix{ComplexF64}(I, nodes, nodes), nodes, nodes, 1, 1)
        actual, roundoff, _ = pi_apply(impulses, energy;
            workers = min(2, Base.Threads.nthreads(:default)))
        coefficients = Dict(lag => coefficient_oracle(lag) for lag in -(nodes-1):(nodes-1))
        exact = [coefficients[i-j] for i in 1:nodes, j in 1:nodes]
        @test reshape(actual, nodes, nodes) ≈ exact rtol=4e-13 atol=2e-14
        @test roundoff < 2e-13
        # Antisymmetry includes the full first and last endpoint hats.
        @test norm(exact+transpose(exact)) == 0
    end
    nodes = 49
    energy = collect(range(-2.0, 3.0; length = nodes))
    plateau = ones(ComplexF64, nodes, 1, 1, 1)
    actual, _, _ = pi_apply(plateau, energy)
    exact = setprecision(256) do
        [Float64((F(BigFloat(i))-F(BigFloat(i-1))-F(BigFloat(nodes-i+1))+
                  F(BigFloat(nodes-i)))/(2big(π))) for i in 1:nodes]
    end
    @test vec(actual) ≈ exact rtol=4e-13 atol=2e-14
    # A represented continuous triangle gives an independent PV at every node.
    energy = collect(range(-4.0, 4.0; length = 33))
    triangle = reshape(ComplexF64[max(1-abs(e), 0) for e in energy], :, 1, 1, 1)
    actual, _, _ = pi_apply(triangle, energy)
    @test vec(actual) ≈ [hat_oracle(e, -1.0, 0.0, 1.0) for e in energy] rtol=4e-13 atol=2e-14
end

@testset "Uniform PI FFT complex dense crosscheck, reuse, geometry and weights" begin
    Random.seed!(0x1335049)
    energy = collect(range(-1.1, 1.4; length = 41))
    Γ = randn(ComplexF64, 41, 3, 2, 2)
    # Include PSD Hermitian blocks with complex off-diagonal entries.
    for e in 1:41, k in 1:3
        z = Matrix(view(Γ, e, k, :, :))
        Γ[e, k, :, :] .= z*z'
    end
    dense = QCLNEGF.product_integration_hilbert_transform(Γ, energy, ones(41))
    for chunk in (1, 3, 32), workers in (1, min(2, Base.Threads.nthreads(:default)))
        actual, roundoff, plan = pi_apply(Γ, energy; chunk, workers)
        @test actual ≈ dense rtol=6e-13 atol=6e-14
        @test isfinite(roundoff) && roundoff < 3e-13
        second = reverse(Γ; dims = 1) .* (0.7+0.2im)
        reused, _ = N._production_pi_fft_hilbert_transform(second, energy, ones(41), plan)
        @test reused ≈ QCLNEGF.product_integration_hilbert_transform(second, energy, ones(41)) rtol=6e-13 atol=6e-14
        reweighted, _ = N._production_pi_fft_hilbert_transform(Γ, energy, 1 .+ abs.(energy), plan)
        @test reweighted == actual
        @test_throws ArgumentError N._production_pi_fft_hilbert_transform(Γ, energy, zeros(41), plan)
    end
    transformed = collect(range(3.7*first(energy)+0.41, 3.7*last(energy)+0.41; length = 41))
    scaled, _, _ = pi_apply(Γ, transformed)
    actual, _, _ = pi_apply(Γ, energy)
    @test scaled ≈ actual rtol=4e-14 atol=4e-14
    # A real geometric perturbation, even one ulp, must never be snapped.
    almost = copy(energy)
    almost[20] = nextfloat(almost[20])
    @test N._uniform_product_integration_geometry(almost) === nothing
    @test_throws ArgumentError N._ProductionPIFFTHilbertPlan(almost, 12)
    # Canonical ranges with inadequately resolved coordinate differences fall back.
    @test N._uniform_product_integration_geometry(collect(range(1e12, 1e12+1; length = 18))) === nothing
end

function production_fixture(nodes; options)
    original = tutorial_numerics()
    changes = (; N_b = 2, N_k = 2, N_E = nodes)
    numerical = NumericalParameters(;
        (name => hasproperty(changes, name) ? getproperty(changes, name) :
                 getfield(original, name) for name in fieldnames(NumericalParameters))...,
    )
    return build_problem_production(;
        numerical, production_options = options,
        scattering = ScatteringOptions(
            LO = false, acoustic = false, impurity = false, IFR = false, alloy = false),
    ).problem
end

function nonuniform_fixture(problem)
    energy = copy(problem.grids.ε)
    energy[25] += (energy[26]-energy[25])/10
    weights = vcat((energy[2]-energy[1])/2, (energy[3:end]-energy[1:end-2])/2,
                   (energy[end]-energy[end-1])/2)
    g = problem.grids
    grids = ModelGrids(g.x, g.wˣ, energy, weights, g.trusted_energy,
                       g.κ, g.wᵏ, g.qᶻ, g.wᑫᶻ, g.φ)
    sp = QCLNEGF._scaled_physics(problem.physical, problem.scales)
    shifts = map(δ -> build_shift_plan(energy, δ),
                 (+sp.Eᵖ, -sp.Eᵖ, +sp.ħωᴸᴼ, -sp.ħωᴸᴼ))
    return NEGFProblem(problem.physical, problem.numerical, problem.scattering,
        problem.scales, grids, problem.profiles, problem.basis, problem.kernels,
        shifts..., :nodal_linear, problem.models)
end

@testset "Selected PI cache memory admission and nonuniform dense fallback" begin
    options = ProductionOptions(worker_count = 1, hilbert_columns = 1,
        algorithms = AlgorithmOptions(hilbert = :product_integration))
    problem = production_fixture(49; options)
    estimate = estimate_production_memory(problem; options)
    cache = build_production_cache(problem; options)
    QCLNEGF._validated_production_cache(cache, problem, options)
    plan_bytes = production_fft_workspace_bytes(cache.product_hilbert_operator) +
                 production_fft_workspace_bytes(cache.product_hilbert_operator; verification = true)
    @test estimate.breakdown[:fft_workspace] == plan_bytes
    # PI always measures arithmetic sensitivity, including when the nodal FFT
    # verification switch is false, as the original dense selected path did.
    unchecked_fft = QCLNEGF.with_production_options(options; verify_fft_roundoff = false)
    @test estimate_production_memory(problem; options = unchecked_fft).peak_bytes == estimate.peak_bytes
    @test estimate.breakdown[:product_hilbert_operator] == 0
    @test estimate.breakdown[:product_hilbert_verification_temporary] == 0
    Γ = randn(ComplexF64, 49, 2, 2, 2)
    retarded, _, witness = QCLNEGF._retarded_self_energy_selected(
        0.3im .* Γ, -0.7im .* Γ, problem, cache, options; return_roundoff = true)
    dense = QCLNEGF.product_integration_hilbert_transform(Γ, problem.grids.ε, problem.grids.wᴱ)
    @test retarded ≈ dense - 0.5im .* Γ rtol=7e-13 atol=7e-14
    @test isfinite(witness) && witness < 3e-13
    @test_throws ArgumentError QCLNEGF._validated_production_cache(cache, problem,
        QCLNEGF.with_production_options(options; hilbert_columns = 2))

    nonuniform = nonuniform_fixture(problem)
    dense_estimate = estimate_production_memory(nonuniform; options)
    @test dense_estimate.breakdown[:fft_workspace] == 0
    @test dense_estimate.breakdown[:product_hilbert_operator] == 49^2*sizeof(Float64)
    @test dense_estimate.breakdown[:product_hilbert_verification_temporary] == 49^2*sizeof(Float64)
    dense_cache = build_production_cache(nonuniform; options)
    impulse = zeros(ComplexF64, 49, 2, 2, 2)
    impulse[25, :, 1, 1] .= 1
    selected, _, _ = QCLNEGF._retarded_self_energy_selected(
        0.3im .* impulse, -0.7im .* impulse, nonuniform, dense_cache, options; return_roundoff = true)
    reference = QCLNEGF.product_integration_hilbert_transform(impulse,
        nonuniform.grids.ε, nonuniform.grids.wᴱ) - 0.5im .* impulse
    @test selected ≈ reference rtol=6e-13 atol=6e-14
    e = nonuniform.grids.ε
    @test real(selected[25, 1, 1, 1]) ≈ log((e[25]-e[24])/(e[26]-e[25]))/(2π) atol=2e-14

    # At S01's energy count, actual cache construction fits a budget smaller
    # than one dense H. No SCBA state or large dense reference is constructed.
    limited = QCLNEGF.with_production_options(options; memory_budget_bytes = 512*1024^2)
    large = production_fixture(30721; options = limited)
    large_cache = build_production_cache(large; options = limited)
    @test large_cache.estimate.peak_bytes < limited.memory_budget_bytes < 30721^2*sizeof(Float64)
    large_plan_bytes = production_fft_workspace_bytes(large_cache.product_hilbert_operator) +
        production_fft_workspace_bytes(large_cache.product_hilbert_operator; verification = true)
    @test large_plan_bytes == 3*65536*2*sizeof(ComplexF64)
    @test large_cache.estimate.breakdown[:fft_workspace] == large_plan_bytes
    too_small = QCLNEGF.with_production_options(options;
        memory_budget_bytes = estimate.peak_bytes-1)
    @test_throws ArgumentError build_production_cache(problem; options = too_small)
    @test_throws ArgumentError QCLNEGF._validated_production_cache(cache, problem, too_small)
end
end
