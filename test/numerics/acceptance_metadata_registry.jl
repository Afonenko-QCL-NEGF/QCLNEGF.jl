module AcceptanceMetadataRegistry
include("../support/common.jl")
const BN=QCLNEGF
const Numerics=BN.QCLNumerics

function assessed_fixture(
    problem;
    algorithms = AlgorithmOptions(),
    context = Dict(),
    metrics = Dict{Symbol,Float64}(),
)
    h=project_hamiltonians(problem, zeros(problem.numerical.N_z))
    green, _=BN._seed_green(problem, h)
    shape=size(green.Gᴿ)
    empty_family=SelfEnergyFamily(
        zeros(ComplexF64, shape),
        zeros(ComplexF64, shape),
        zeros(ComplexF64, shape),
    )
    scba=SCBAResult(
        green,
        Dict{Symbol,SelfEnergyFamily}(),
        empty_family,
        empty_family,
        empty_family,
        SCBAIteration[],
        false,
        :max_iterations,
        :unresolved,
    )
    options=SolverOptions()
    observables=Dict{Symbol,Any}(
        :restart_contract=>BN._solver_restart_contract(options, algorithms),
        :acceptance_measurement_context=>context,
    )
    return NEGFSolution(
        problem,
        options,
        zeros(problem.numerical.N_z),
        zeros(problem.numerical.N_z),
        scba,
        OuterIteration[],
        observables,
        ConvergenceReport(false, metrics, String[]),
        false,
        :max_iterations,
    )
end

@testset "Every final acceptance metric owns a complete descriptive contract" begin
    problem=build_problem(
        numerical = tutorial_numerics(),
        physical = reference_parameters(Tᴸ = 70u"K"),
        scattering = ScatteringOptions(
            LO = false,
            acoustic = false,
            impurity = false,
            IFR = false,
            alloy = false,
        ),
    )
    solution=assessed_fixture(problem)
    extra=Dict{Symbol,Float64}(
        name=>0.0 for name in (
            :r_energy_electron_electron,
            :r_hot_lo_power,
            :r_D_candidate,
            :r_A_candidate,
            :r_lambda_candidate,
        )
    )
    for mechanism in (:LO, :acoustic, :impurity, :IFR, :alloy, :electron_electron),
        prefix in (:r_C_, :r_ImC_, :r_roundC_)

        extra[Symbol(prefix, mechanism)]=0.0
    end
    limits=BN._solution_validation_limits(solution.options.tolerances, extra)
    for (name, threshold) in limits
        metadata=Numerics._acceptance_metric_metadata(name, solution)
        @test metadata["registry_version"]=="qcl-negf-acceptance-metadata-v1"
        for key in (
            "formula",
            "formula_source",
            "formula_context",
            "units",
            "normalization",
            "weights",
        )
            @test metadata[key] isa String && !isempty(metadata[key])
        end
        @test Set(keys(metadata["denominator_floor"]))==Set(("kind", "value", "units"))
        @test metadata["threshold_origin"]["kind"]=="project_acceptance_policy"
        @test !metadata["threshold_origin"]["literature_threshold"]
        @test metadata["error_budget"]["bound_status"]=="not_measured"
        @test metadata["error_budget"]["estimated_bound"]===nothing
        @test metadata["applicability"]["applies"] isa Bool
        @test !isempty(metadata["applicability"]["reason"])
        # Every machine reference resolves to a real local evaluator file.
        for reference in split(metadata["formula_source"], ';')
            path=strip(first(split(reference, "::")))
            @test isfile(joinpath(TEST_ROOT, "..", path))
        end
    end
    @test_throws ArgumentError Numerics._acceptance_metric_metadata(
        :undefined_future_metric,
        solution,
    )
    assessment=solution_scientific_assessment(solution)
    @test assessment["metrics"]["r_K"]["status"]=="not_measured"
    @test assessment["metrics"]["r_edge_gamma"]["status"]=="not_applicable"
    @test assessment["metrics"]["r_roundoff"]["status"]=="not_applicable"
    @test !assessment["scientific_accepted"]
    @test assessment["quadrature"]["energy"]["weights"]==problem.grids.wᴱ
    @test assessment["scales"]["E0_eV"]==problem.scales.E₀_eV

    # Match the actual final evaluator at the equality boundary. Iteration
    # confirmation still has a separate strict-inequality contract.
    equal=assessed_fixture(problem; metrics = Dict(:r_sum=>SolverTolerances().r_sum))
    record=solution_scientific_assessment(equal)["metrics"]["r_sum"]
    @test record["status"]=="pass"
    @test record["comparison"]=="less_or_equal"
    @test record["denominator_floor"]["kind"]=="none"
    bad=assessed_fixture(problem; metrics = Dict(:r_sum=>NaN))
    @test solution_scientific_assessment(bad)["metrics"]["r_sum"]["status"]=="error"

    # Product integration really uses max(norm,eps), unlike the additive
    # 1e-14 of direct/FFT. Metadata must describe the selected operator.
    product=assessed_fixture(
        problem;
        algorithms = AlgorithmOptions(hilbert = :product_integration),
    )
    product_metadata=Numerics._acceptance_metric_metadata(:r_roundoff, product)
    @test product_metadata["denominator_floor"]["kind"]=="maximum"
    @test product_metadata["denominator_floor"]["value"]==eps(Float64)
    @test Numerics._acceptance_metric_metadata(:r_PSD, solution)["denominator_floor"]["value"]==floatmin(
        Float64,
    )
    @test Numerics._acceptance_metric_metadata(:r_Jchange, solution)["denominator_floor"]["units"]=="A/m^2"
    @test Numerics._acceptance_metric_metadata(:r_J, solution)["denominator_floor"]["units"]=="scaled boundary number flux"

    active=build_problem(
        numerical = tutorial_numerics(),
        physical = reference_parameters(Tᴸ = 70u"K"),
        scattering = ScatteringOptions(
            LO = false,
            acoustic = true,
            impurity = false,
            IFR = false,
            alloy = false,
        ),
    )
    unchecked=assessed_fixture(
        active;
        context = Dict("hilbert_roundoff_measured"=>false),
        metrics = Dict(:r_roundoff=>0.0),
    )
    unchecked_record=solution_scientific_assessment(unchecked)["metrics"]["r_roundoff"]
    @test unchecked_record["status"]=="not_measured"
    @test unchecked_record["value"]===nothing
    @test unchecked_record["recorded_value"]==0.0
    @test !solution_scientific_assessment(unchecked)["physical_gates_passed"]
end
end
