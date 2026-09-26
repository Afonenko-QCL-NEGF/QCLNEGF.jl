module ProductionBoundedOperations
include("../support/common.jl")
include("../support/native_physics_fixture.jl")
const BN = QCLNEGF
const N = BN.QCLNumerics

@testset "Bounded cavity agrees with independent full finite-chain oracle" begin
    fixture = native_physics_fixture(; energy_nodes = 17)
    problem, h = fixture.problem, project_hamiltonians(fixture.problem, fixture.Uᴴ)
    ne, nk, nb = problem.numerical.N_E, problem.numerical.N_k, problem.numerical.N_b
    family = SelfEnergyFamily(
        zeros(ComplexF64, ne, nk, nb, nb),
        zeros(ComplexF64, ne, nk, nb, nb),
        zeros(ComplexF64, ne, nk, nb, nb),
    )
    for e = 1:ne, m = 1:nk, a = 1:nb
        width = 0.02 + 0.001e + 0.002m
        family.Σᴿ[e, m, a, a] = -im*width
        family.Σˡ[e, m, a, a] = 0.6im*width
        family.Σᵍ[e, m, a, a] = -1.4im*width
    end
    for shift in (:dense, :sparse_plan, :conservative_pair),
        periods in (1, 2),
        workers in (1, 2)

        options = ProductionOptions(
            worker_count = workers,
            algorithms = AlgorithmOptions(
                embedding = :finite_chain,
                embedding_periods = periods,
                energy_shift = shift,
            ),
        )
        expected = BN.cavity_embedding_self_energy(
            problem,
            family,
            h;
            periods,
            energy_shift = shift,
            worker_count = 1,
        )
        actual = N._cavity_embedding_production(problem, family, h, options)
        for index in eachindex(expected), component in (:Σᴿ, :Σˡ, :Σᵍ)
            @test getfield(actual[index], component) ≈ getfield(expected[index], component) rtol=5e-12 atol=2e-13
        end
    end
end

@testset "Checkpoint decisions precede snapshots and preserve marker cadence" begin
    fixture = native_physics_fixture(; energy_nodes = 17)
    calls, requests = Ref(0), Int[]
    options = SolverOptions(
        max_scba = 3,
        max_poisson = 1,
        convergence = ConvergencePolicy(
            required_consecutive_scba_passes = 99,
            stagnation_window = 0,
            stagnation_relative_improvement = 0.0,
        ),
    )
    production = ProductionOptions(
        worker_count = 1,
        checkpoint_every_scba = 1,
        checkpoint_request = context -> (push!(requests, context.iteration); false),
        physics_markers = SCBAPhysicsMarkerPolicy(cadence = 10),
        progress_every_scba = 0,
    )
    result = solve_scba_production(
        fixture.problem,
        fixture.Uᴴ;
        options,
        production_options = production,
        iteration_callback = state->(calls[]+=1),
    )
    @test requests == [1, 2, 3]
    @test calls[] == 0
    @test result.history[2].physical_markers === nothing
    @test result.history[end].physical_markers !== nothing
end

@testset "Producer phases negotiate bounded task widths and actual counter intervals" begin
    events = N.SolverEvent[]
    options = ProductionOptions(
        phase_timing = true,
        event_sink = e->push!(events, e),
        phase_request = context->min(2, context.task_width),
    )
    phase = N._production_phase_begin(options, :dyson, 3; task_width = 4, work_units = 12)
    @test phase.options.worker_count == 2
    scratch = ones(1000)
    @test sum(scratch)==1000
    N._production_phase_end(options, phase)
    @test [e.action for e in events] == [:phase_begin, :phase_end]
    values = Dict(m.name=>m.value for m in events[end].metrics)
    @test values[:end_monotonic_ns] >= values[:start_monotonic_ns]
    @test values[:allocated_bytes] >= 0
    @test values[:gc_seconds] >= 0
    @test values[:task_width] == 2
    @test_throws ArgumentError N._production_phase_begin(
        options,
        :residuals,
        3;
        task_width = 4,
        elastic = false,
    )
end

@testset "Anderson history recycles buffers and caches the exact Gram matrix" begin
    family(value) = SelfEnergyFamily(
        fill(ComplexF64(value), 2, 1, 1, 1),
        zeros(ComplexF64, 2, 1, 1, 1),
        zeros(ComplexF64, 2, 1, 1, 1),
    )
    workspace = N._AndersonWorkspace()
    ids = Set{UInt}()
    for step = 1:8
        N._anderson_append!(workspace, [family(0)], [family(1/step)], 3)
        for state in workspace.states
            push!(ids, objectid(only(state).Σᴿ))
        end
        @test length(workspace.states) <= 3
        for i in eachindex(workspace.residuals), j in eachindex(workspace.residuals)
            @test workspace.gram[i, j] ==
                  N._anderson_inner(workspace.residuals[i], workspace.residuals[j])
        end
    end
    @test length(ids)==3
    # A reset reuses the same bounded allocation pool on subsequent steps.
    N._anderson_append!(workspace, [family(0)], [family(2)], 3)
    @test length(workspace.states)==1
    @test length(workspace.states)+length(workspace.free_states)==3
    @test N._anderson_guard([family(1)], ProductionOptions(worker_count = 2)).passed
end
@testset "A cached FFT plan obeys a smaller task allocation without changing arithmetic" begin
    energy = collect(range(-1.0, 1.0; length = 17))
    weights = fill(0.125, 17)
    weights[[1, end]] ./= 2
    gamma = reshape(ComplexF64[sin(i/10) for i = 1:(17*2*2*2)], 17, 2, 2, 2)
    plan = ProductionFFTHilbertPlan(energy, 8; worker_count = 2, verify_roundoff = true)
    serial, serial_error = production_fft_hilbert_transform(
        gamma,
        energy,
        weights,
        plan;
        return_roundoff = true,
        worker_count = 1,
    )
    parallel, parallel_error = production_fft_hilbert_transform(
        gamma,
        energy,
        weights,
        plan;
        return_roundoff = true,
        worker_count = 2,
    )
    @test serial == parallel
    @test serial_error == parallel_error
end

@testset "Parallel Anderson safeguard preserves the literal matrix tests" begin
    function literal(families)
        for f in families, e in axes(f.Σᴿ, 1), m in axes(f.Σᴿ, 2)
            r, l, g = (Matrix(view(getfield(f, c), e, m, :, :)) for c in (:Σᴿ, :Σˡ, :Σᵍ))
            all(isfinite, r) && all(isfinite, l) && all(isfinite, g) || return false
            eigmax(Hermitian((r-r')/(2im))) <=
            64eps(Float64)*max(opnorm(r), floatmin(Float64)) || return false
            for block in (-im .* l, im .* g)
                budget = psd_error_budget(opnorm(block), size(block, 1))
                norm(block-block')/2 <= budget || return false
                eigmin(Hermitian((block+block')/2)) >= -budget || return false
            end
            norm(r-r'-(g-l)) <=
            psd_error_budget(max(opnorm(r), opnorm(l), opnorm(g)), size(r, 1)) ||
                return false
        end
        true
    end
    family = SelfEnergyFamily((zeros(ComplexF64, 8, 2, 3, 3) for _ = 1:3)...)
    for e = 1:8, m = 1:2
        b = ComplexF64[sin(e+a*j)+im*cos(m+a+2j) for a = 1:3, j = 1:3]
        c = ComplexF64[cos(e+a*j)+im*sin(m+2a+j) for a = 1:3, j = 1:3]
        incoming, outgoing = b*b', c*c'
        h = ComplexF64[sin(a+j)+im*(a-j)/3 for a = 1:3, j = 1:3]
        family.Σᴿ[e, m, :, :] .= h - 0.5im .* (incoming+outgoing)
        family.Σˡ[e, m, :, :] .= im .* incoming
        family.Σᵍ[e, m, :, :] .= -im .* outgoing
    end
    for defect in (0.0, 1e-6, -1e-6), workers in (1, 2)
        altered = BN._copy_family(family)
        altered.Σˡ[3, 2, 1, 1] += defect
        @test N._anderson_guard([altered], ProductionOptions(worker_count = workers)).passed ==
              literal([altered])
    end
end

end
