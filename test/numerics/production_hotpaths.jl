module ProductionHotpaths
include("../support/common.jl")
include("../support/production_residuals.jl")
include("../support/native_physics_fixture.jl")
const N = QCLNEGF.QCLNumerics

@testset "Tiny inverse retains pivoting, scaling and singular detection" begin
    Random.seed!(1082)
    for n = 1:9
        workspace = N._ProductionInverseWorkspace(n)
        random = randn(ComplexF64, n, n)
        q = Matrix(qr(random).Q)
        # Non-normal and near-resonant retarded matrices exercise all pivots.
        cases = [
            random + 0.1I,
            q*Diagonal(n == 1 ? [1e-10] : collect(range(1e-10, 1; length = n)))*q' +
            1e-12im*I,
        ]
        if n > 1
            pivoted = Matrix{ComplexF64}(I, n, n)
            pivoted[[1, n], :] = pivoted[[n, 1], :]
            push!(cases, pivoted)
        end
        for input in cases, scale in (1e-150, 1.0, 1e150)
            a = input*scale
            output = N._cavity_inverse!(copy(a), workspace)
            reference = inv(a)
            # Backward error is the appropriate stability measure near resonance.
            denominator = norm(a)*norm(output) + sqrt(n)
            @test norm(a*output-I)/denominator < 32n*eps(Float64)
            @test norm(output-reference)/norm(reference) < 64n*eps(Float64)*cond(a)
            if cond(a) < 1e4
                @test output ≈ reference rtol=2e-13 atol=0
            end
        end
        @test_throws SingularException N._cavity_inverse!(
            zeros(ComplexF64, n, n),
            workspace,
        )
        if n > 1
            deficient = Matrix{ComplexF64}(I, n, n)
            deficient[end, :] .= deficient[1, :]
            @test_throws SingularException N._cavity_inverse!(deficient, workspace)
        end
        @test_throws DimensionMismatch N._cavity_inverse!(
            zeros(ComplexF64, n, n+1),
            workspace,
        )
    end
end

@testset "Leased residual width preserves every exact reduction" begin
    ε, h, green, total, candidate = _random_residual_fixture(; NE = 17, Nk = 3, Nb = 5)
    workspace = ProductionResidualWorkspace(17, 3, 5; chunk_size = 3)
    expected = production_residual_suite(
        ε,
        h,
        green,
        total,
        candidate;
        workspace,
        worker_count = 1,
    )
    expected_sigma =
        production_selfenergy_residual(total, candidate; workspace, worker_count = 1)
    for workers in (0, 1, 2, 31)
        @test production_residual_suite(
            ε,
            h,
            green,
            total,
            candidate;
            workspace,
            worker_count = workers,
        ) == expected
        @test production_selfenergy_residual(
            total,
            candidate;
            workspace,
            worker_count = workers,
        ) == expected_sigma
    end
    @test N._production_residual_worker_count(workspace, 17, 1) == 1
    @test N._production_residual_worker_count(workspace, 17, 31) ==
          min(31, Threads.nthreads(:default))
    @test_throws ArgumentError production_residual_suite(
        ε,
        h,
        green,
        total,
        candidate;
        workspace,
        worker_count = -1,
    )
    @test_throws ArgumentError production_selfenergy_residual(
        total,
        candidate;
        workspace,
        worker_count = -1,
    )
    causal = SelfEnergyFamily(copy(total.Σᴿ), copy(total.Σˡ), copy(total.Σᵍ))
    for e = 1:17, m = 1:3
        block = randn(ComplexF64, 5, 5)
        causal.Σᴿ[e, m, :, :] .= block+block' - im*(block*block'+I)
    end
    result = production_residual_suite(ε, h, green, causal, candidate; workspace)
    @test result.r_caus == QCLNEGF._causality_residual(causal) == 0.0
end

@testset "Five-state cavity matches the independent finite-chain oracle" begin
    base = tutorial_numerics()
    changes = (N_z = 65, N_b = 5, N_E = 17, N_k = 3, N_qz = 9, N_φ = 8)
    numerical = NumericalParameters(;
        (
            name => (
                hasproperty(changes, name) ? getproperty(changes, name) :
                getfield(base, name)
            ) for name in fieldnames(NumericalParameters)
        )...,
    )
    problem = build_problem(;
        numerical,
        validate_static = false,
        scattering = ScatteringOptions(
            LO = false,
            acoustic = false,
            impurity = false,
            IFR = false,
            alloy = false,
        ),
    )
    h = project_hamiltonians(problem, zeros(numerical.N_z))
    family = SelfEnergyFamily((zeros(ComplexF64, 17, 3, 5, 5) for _ = 1:3)...)
    b = randn(ComplexF64, 5, 5)
    gamma = 0.04*(b*b'+I)
    for e = 1:17, m = 1:3
        family.Σᴿ[e, m, :, :] .= -0.5im*gamma
        family.Σˡ[e, m, :, :] .= 0.3im*gamma
        family.Σᵍ[e, m, :, :] .= -0.7im*gamma
    end
    for shift in (:dense, :sparse_plan, :conservative_pair), periods in (1, 4)
        oracle = QCLNEGF.cavity_embedding_self_energy(
            problem,
            family,
            h;
            periods,
            energy_shift = shift,
            worker_count = 1,
        )
        for workers in (1, 2)
            options = ProductionOptions(
                worker_count = workers,
                algorithms = AlgorithmOptions(
                    embedding = :finite_chain,
                    embedding_periods = periods,
                    energy_shift = shift,
                ),
            )
            actual = N._cavity_embedding_production(problem, family, h, options)
            for (left, right) in zip(actual, oracle), key in (:Σᴿ, :Σˡ, :Σᵍ)
                @test getfield(left, key) ≈ getfield(right, key) rtol=5e-12 atol=2e-13
            end
        end
    end
end

@testset "SCBA negotiates residual width and consumes only explicit private restart" begin
    fixture = native_physics_fixture(; energy_nodes = 17)
    options = SolverOptions(
        max_scba = 3,
        max_poisson = 1,
        convergence = ConvergencePolicy(
            required_consecutive_scba_passes = 99,
            stagnation_window = 0,
            stagnation_relative_improvement = 0.0,
        ),
    )
    events = N.SolverEvent[]
    requests = []
    production = ProductionOptions(
        worker_count = 0,
        checkpoint_every_scba = 1,
        progress_every_scba = 0,
        phase_timing = true,
        event_sink = event->push!(events, event),
        phase_request = context->(push!(requests, context); 1),
    )
    saved = Ref{Any}(nothing)
    uninterrupted = solve_scba_production(
        fixture.problem,
        fixture.Uᴴ;
        options,
        production_options = production,
        iteration_callback = state->(
            last(state.history).ν == 1 && (saved[]=deepcopy(state))
        ),
    )
    @test saved[] !== nothing
    @test all(request.elastic for request in requests if request.name === :residuals)
    widths = [
        only(m.value for m in event.metrics if m.name===:task_width) for
        event in events if event.action===:phase_begin && event.label=="residuals"
    ]
    @test !isempty(widths) && all(==(1), widths)
    original = deepcopy(saved[])
    resumed = solve_scba_production(
        fixture.problem,
        fixture.Uᴴ;
        options,
        production_options = production,
        initial = saved[],
        resume_iterations = true,
    )
    @test saved[].embedding.Σᴿ == original.embedding.Σᴿ
    @test saved[].embedding_plus.Σᴿ == original.embedding_plus.Σᴿ
    private = deepcopy(saved[])
    original_embedding = private.embedding.Σᴿ
    empty!(requests)
    consumed = solve_scba_production(
        fixture.problem,
        fixture.Uᴴ;
        options,
        production_options = production,
        initial = private,
        resume_iterations = true,
        consume_initial = true,
    )
    @test consumed.embedding.Σᴿ === original_embedding
    @test consumed.green.Gᴿ ≈ resumed.green.Gᴿ rtol=2e-13
    @test consumed.green.Gˡ ≈ uninterrupted.green.Gˡ rtol=2e-13
    @test [row.r_Σ for row in consumed.history] ≈ [row.r_Σ for row in uninterrupted.history] rtol=2e-13
    @test first(requests).memory_burst_bytes == 0
    @test :restart_candidate ∉ [r.name for r in requests]
    @test all(
        name in [r.name for r in requests] for
        name in (:restart_scattering, :restart_embedding, :restart_mixing)
    )
    wrapper_private = deepcopy(saved[])
    wrapper_embedding = wrapper_private.embedding.Σᴿ
    wrapper = solve_adaptive_production(
        fixture.problem;
        options,
        production_options = production,
        initial_Uᴴ = fixture.Uᴴ,
        initial_scba = wrapper_private,
        resume_scba = true,
        consume_initial = true,
    )
    @test wrapper.scba.embedding.Σᴿ === wrapper_embedding
    @test wrapper.scba.green.Gᴿ ≈ consumed.green.Gᴿ rtol=2e-13
end
end
