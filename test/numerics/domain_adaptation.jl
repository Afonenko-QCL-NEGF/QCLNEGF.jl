module DomainAdaptationContracts
include("../support/common.jl")

@testset "Declared energy expansion rebuilds operators with a bounded cold seed" begin
    @test_throws ArgumentError DomainAdaptationPolicy(mode = :shrink)
    @test_throws ArgumentError DomainAdaptationPolicy(maximum_energy_nodes = 2)
    base=tutorial_numerics()
    changes=(
        N_z = 25,
        N_b = 2,
        N_E = 17,
        N_k = 2,
        N_qz = 9,
        N_φ = 8,
        E_min = -0.1u"eV",
        E_max = 0.05u"eV",
    )
    n=NumericalParameters(;
        (
            field=>(
                hasproperty(changes, field) ? getproperty(changes, field) :
                getfield(base, field)
            ) for field in fieldnames(NumericalParameters)
        )...,
    )
    algorithms=AlgorithmOptions(
        energy_shift = :dense,
        hilbert = :direct,
        contraction = :dense_blas,
        kernel_build = :direct,
        localization = :none,
    )
    # Test domain rebuilding on an equilibrium problem with a physical bath.
    # Its LO linewidth gives the sampled Keldysh pair finite spectral capacity.
    p=reference_parameters(Tᴸ = 70u"K", V_period = 0u"mV")
    scattering=ScatteringOptions(
        LO = true,
        acoustic = false,
        impurity = false,
        IFR = false,
        alloy = false,
    )
    kernel_options=ProductionKernelOptions()
    source=build_configured_scattering_problem(
        physical = p,
        numerical = n,
        scales = ScaleSystem(),
        scattering = scattering,
        algorithms = algorithms,
        kernel_options = kernel_options,
    ).problem
    original=representation_coverage(source)
    @test !original.covered
    production=ProductionOptions(
        algorithms = algorithms,
        parallel_backend = :blas,
        worker_count = 1,
        progress_every_scba = 0,
        progress_every_outer = 0,
        checkpoint_every_scba = 0,
        checkpoint_every_outer = 0,
    )
    policy=DomainAdaptationPolicy(
        mode = :expand_energy_window,
        maximum_expansions = 1,
        maximum_energy_nodes = 10_000,
        tail_threshold = 1e-6,
    )
    result=solve_adaptive_production(
        source;
        domain_adaptation = policy,
        kernel_options,
        options = SolverOptions(max_scba = 2, max_poisson = 1),
        production_options = production,
    )
    @test result.problem.models == source.models
    @test result.problem.grids !== source.grids
    @test result.problem.kernels !== source.kernels
    @test representation_coverage(result.problem).discretization_id !=
          original.discretization_id
    @test result.problem.numerical.N_E > source.numerical.N_E
    @test maximum(diff(result.problem.grids.ε)) <= maximum(diff(source.grids.ε))*(1+1e-12)
    record=result.observables[:domain_adaptation]
    @test record["expansions"] == 1
    @test length(record["revisions"]) == 2
    @test record["revisions"][1]["source"] == "preflight_estimate"
    @test record["revisions"][2]["source"] == "measured_state"
    @test !result.converged
    @test_throws ArgumentError solve_adaptive_production(
        source;
        domain_adaptation = DomainAdaptationPolicy(
            mode = :expand_energy_window,
            maximum_energy_nodes = 3,
        ),
        kernel_options,
        options = SolverOptions(max_scba = 1, max_poisson = 1),
        production_options = production,
    )
end
end
