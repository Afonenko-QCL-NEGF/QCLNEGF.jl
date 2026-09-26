module ProductionExactOperations
include("../support/common.jl")
include("../support/native_physics_fixture.jl")
const BN = QCLNEGF
const N = BN.QCLNumerics

@testset "Owned LU inverse matches the dense oracle and rejects singular blocks" begin
    for nb in (1, 2, 5, 8, 16)
        matrix = ComplexF64[sin(a+2b)+im*cos(2a-b) for a = 1:nb, b = 1:nb] + 4I
        work = N._ProductionInverseWorkspace(nb)
        output = copy(matrix)
        @test N._production_inverse!(output, work) === output
        @test output ≈ inv(matrix) rtol=5e-13 atol=5e-14
        @test matrix * output ≈ Matrix{ComplexF64}(I, nb, nb) rtol=5e-13 atol=5e-13
        copyto!(output, matrix)
        N._production_inverse!(output, work) # warm this exact specialization
        copyto!(output, matrix)
        @test @allocated(N._production_inverse!(output, work)) <= 1024
    end
    @test_throws SingularException N._production_inverse!(
        zeros(ComplexF64, 2, 2),
        N._ProductionInverseWorkspace(2),
    )
    @test_throws DimensionMismatch N._production_inverse!(
        ones(ComplexF64, 3, 3),
        N._ProductionInverseWorkspace(2),
    )
end

@testset "Production observables preserve independent matrix-product definitions" begin
    fixture = native_physics_fixture(; energy_nodes = 17)
    problem = fixture.problem
    shape = size(fixture.scba.green.Gᴿ)
    field(phase) =
        reshape(ComplexF64[sin(i/7+phase)+im*cos(i/11-phase) for i = 1:prod(shape)], shape)
    gr, gl, gg, a = (field(phase) for phase in (0.1, 0.3, 0.7, 1.1))
    green = BN.GreenState(gr, gl, gg, a, ones(shape[1:2]), ones(shape[1:2]))
    family = SelfEnergyFamily(field(1.3), field(1.7), field(1.9))
    for e in (1, 8, 17), m in (1, 3)
        expected = BN._current_trace(family.Σᵍ, family.Σˡ, gl, gg, e, m)
        actual = N._production_current_trace(family.Σᵍ, family.Σˡ, gl, gg, e, m)
        @test actual ≈ expected rtol=2e-13 atol=2e-13
    end
    @test N._production_boundary_flux_bar(problem, family, green) ≈
          BN._boundary_flux_bar(problem, family, green) rtol=3e-13 atol=2e-13
    target = BN._scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    populations = real.(diag(BN._sheet_density_matrix_bar(problem, gl))) ./ target
    @test N._production_state_populations(problem, gl) ≈ populations rtol=3e-13 atol=2e-13
    measurement =
        N._production_observable_measurement(problem, green, family, nothing, nothing)
    @test measurement.populations ≈ populations rtol=3e-13 atol=2e-13
    @test measurement.rJchange == Inf
    @test measurement.rpopulation == Inf
end

@testset "Reusable shift output preserves each discretization and rejects shared storage" begin
    fixture = native_physics_fixture(; energy_nodes = 17)
    for shift in (:dense, :sparse_plan, :conservative_pair)
        problem = retarget_problem(fixture.problem; energy_shift = shift)
        options = ProductionOptions(
            worker_count = 1,
            verify_fft_roundoff = false,
            algorithms = AlgorithmOptions(energy_shift = shift),
        )
        cache = build_production_cache(problem; options)
        input = copy(fixture.scba.green.Gˡ)
        before = copy(input)
        output = similar(input)
        for direction in (:plus_lo, :minus_lo)
            matrix = direction === :plus_lo ? problem.W₊ᴸᴼ : problem.W₋ᴸᴼ
            expected = apply_energy_shift(matrix, input)
            @test N._selected_energy_shift!(
                output,
                problem,
                cache,
                direction,
                input,
                options,
            ) === output
            @test output ≈ expected rtol=5e-14 atol=5e-14
            @test input == before
        end
        @test_throws ArgumentError N._selected_energy_shift!(
            input,
            problem,
            cache,
            :plus_lo,
            input,
            options,
        )
        alias = reshape(view(vec(input), :), size(input))
        @test_throws ArgumentError N._selected_energy_shift!(
            alias,
            problem,
            cache,
            :plus_lo,
            input,
            options,
        )
        @test_throws ArgumentError N._selected_energy_shift!(
            output,
            problem,
            cache,
            :invalid,
            input,
            options,
        )
    end
end

@testset "SCBA-scoped cavity plan is reusable only for its declared Hartree problem" begin
    fixture = native_physics_fixture(; energy_nodes = 17)
    problem = fixture.problem
    h = project_hamiltonians(problem, fixture.Uᴴ)
    shape = size(fixture.scba.green.Gᴿ)
    family = SelfEnergyFamily((zeros(ComplexF64, shape) for _ = 1:3)...)
    for e in axes(family.Σᴿ, 1), m in axes(family.Σᴿ, 2), a in axes(family.Σᴿ, 3)
        family.Σᴿ[e, m, a, a] = -0.02im
        family.Σˡ[e, m, a, a] = 0.01im
        family.Σᵍ[e, m, a, a] = -0.03im
    end
    for workers in (1, 2)
        options = ProductionOptions(
            worker_count = workers,
            algorithms = AlgorithmOptions(embedding = :finite_chain, embedding_periods = 2),
        )
        plan = N._cavity_production_plan(problem, h, options)
        expected = BN.cavity_embedding_self_energy(
            problem,
            family,
            h;
            periods = 2,
            energy_shift = options.algorithms.energy_shift,
            worker_count = 1,
        )
        for repeat_index = 1:2
            actual = N._cavity_embedding_production(problem, family, h, options; plan)
            for index in eachindex(actual), component in (:Σᴿ, :Σˡ, :Σᵍ)
                @test getfield(actual[index], component) ≈
                      getfield(expected[index], component) rtol=5e-12 atol=2e-13
            end
        end
        altered = copy(h)
        altered[1, 1, 1] += 1e-4
        @test_throws ArgumentError N._cavity_embedding_production(
            problem,
            family,
            altered,
            options;
            plan,
        )
    end
end

@testset "Memory planning uses explicit launch capacity independently of planner threads" begin
    numerical = tutorial_numerics()
    options = ProductionOptions(
        worker_count = 1,
        energy_chunk = 8,
        algorithms = AlgorithmOptions(embedding = :finite_chain),
    )
    one = estimate_production_memory(numerical; options, worker_capacity = 1)
    eight = estimate_production_memory(numerical; options, worker_capacity = 8)
    @test eight.peak_bytes >= one.peak_bytes
    @test eight.breakdown[:blas_workspace] > one.breakdown[:blas_workspace]
    @test eight.breakdown[:cavity_workspace] > one.breakdown[:cavity_workspace]
    @test_throws ArgumentError estimate_production_memory(
        numerical;
        options,
        worker_capacity = 0,
    )
end

end
