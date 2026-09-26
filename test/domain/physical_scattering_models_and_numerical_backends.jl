module Suite_T027
include("../support/common.jl")
include("../support/configuration_schema_scattering_domain.jl")

@testset "Physical scattering models and numerical backends" begin
    physical = reference_parameters()
    models = QCLNEGF.scattering_models(physical, default_scattering())
    @test QCLNEGF.scattering_id.(models) == [:LO, :acoustic, :impurity, :IFR]
    @test QCLNEGF.scattering_options(models) == default_scattering()
    @test QCLNEGF._models_match_physical(models, physical)
    @test all(model -> QCLNEGF.scattering_model_role(model) === :physical_model, models)

    invalid_lo = QCLNEGF.LOPhononModel(
        0.0u"eV",
        physical.Tᴸᴼ,
        physical.ε_s,
        physical.ε_∞,
        physical.qᴸᴼ_s,
    )
    @test_throws QCLNEGF.ScatteringValidationError begin
        QCLNEGF.validate_scattering_model(invalid_lo)
    end
    @test_throws QCLNEGF.ScatteringValidationError begin
        QCLNEGF.LOPhononModel(
            phonon_energy = 1.0u"m",
            temperature = physical.Tᴸᴼ,
            static_relative_permittivity = physical.ε_s,
            high_frequency_relative_permittivity = physical.ε_∞,
            screening_wavenumber = physical.qᴸᴼ_s,
        )
    end

    literal = QCLNEGF.ScatteringNumericalPlan(QCLNEGF.LiteralScatteringKernelBackend())
    exact_blas = QCLNEGF.ScatteringNumericalPlan(
        QCLNEGF.LiteralScatteringKernelBackend();
        contraction = :dense_blas,
    )
    controlled = QCLNEGF.ScatteringNumericalPlan(
        QCLNEGF.TabulatedScatteringKernelBackend();
        contraction = :dense_blas,
    )
    low_rank = QCLNEGF.ScatteringNumericalPlan(
        QCLNEGF.LiteralScatteringKernelBackend();
        contraction = :low_rank,
        low_rank_relative_tolerance = 1e-6,
        low_rank_maximum_rank = 32,
    )
    @test QCLNEGF.scattering_evidence_class(literal) === :E0
    @test QCLNEGF.scattering_evidence_class(exact_blas) === :E1
    @test QCLNEGF.scattering_evidence_class(controlled) === :E2
    @test QCLNEGF.scattering_evidence_class(low_rank) === :E2
    @test QCLNEGF.scattering_evidence_class(
        AlgorithmOptions(contraction = :dense_blas, kernel_build = :direct),
    ) === :E1
    @test QCLNEGF.scattering_evidence_class(
        AlgorithmOptions(kernel_build = :tabulated),
    ) === :E2
    @test QCLNEGF.scattering_evidence_class(
        AlgorithmOptions(retarded_real_part = :drop),
    ) === :E3

    production_direct_options = AlgorithmOptions(
        solver_backend = :production,
        kernel_build = :direct,
        contraction = :literal,
    )
    production_direct = QCLNEGF.scattering_numerical_plan(production_direct_options)
    @test production_direct.kernel_backend isa
          QCLNEGF.ExactParallelScatteringKernelBackend
    @test QCLNEGF.scattering_evidence_class(production_direct) === :E1

    educational_direct_options = AlgorithmOptions(
        solver_backend = :educational,
        kernel_build = :direct,
        contraction = :literal,
    )
    educational_direct = QCLNEGF.scattering_numerical_plan(educational_direct_options)
    @test educational_direct.kernel_backend isa QCLNEGF.LiteralScatteringKernelBackend
    @test QCLNEGF.scattering_evidence_class(educational_direct) === :E0
    @test_throws QCLNEGF.ScatteringValidationError begin
        QCLNEGF.scattering_numerical_plan(
            production_direct_options;
            solver_backend = :educational,
        )
    end
    @test QCLNEGF.scattering_selection_evidence_class(models, models[1:(end-1)]) === :E3
    @test QCLNEGF.scattering_selection_evidence_class(models, reverse(models)) === :E0
    changed_models = copy(models)
    changed_models[1] = QCLNEGF.LOPhononModel(
        phonon_energy = 37.0u"meV",
        temperature = physical.Tᴸᴼ,
        static_relative_permittivity = physical.ε_s,
        high_frequency_relative_permittivity = physical.ε_∞,
        screening_wavenumber = physical.qᴸᴼ_s,
    )
    @test QCLNEGF.scattering_selection_evidence_class(models, changed_models) === :E3
    @test !QCLNEGF._models_match_physical(changed_models, physical)

    @test_throws QCLNEGF.ScatteringValidationError begin
        QCLNEGF.ScatteringNumericalPlan(
            QCLNEGF.LiteralScatteringKernelBackend();
            contraction = :dense_blas,
            low_rank_relative_tolerance = 1e-4,
        )
    end

    numerical = tutorial_numerics()
    shape = (
        numerical.N_k,
        numerical.N_k,
        numerical.N_b,
        numerical.N_b,
        numerical.N_b,
        numerical.N_b,
    )
    normalized = fill(1.0 + 0.0im, shape)
    form_factor = zeros(ComplexF64, numerical.N_qz, numerical.N_b, numerical.N_b)
    acoustic = only(filter(model -> QCLNEGF.scattering_id(model) === :acoustic, models))
    kernels = KernelSet(
        Dict(:acoustic => normalized),
        Dict(:acoustic => 2.0),
        form_factor,
        [:acoustic],
    )
    backend_result = QCLNEGF.ScatteringKernelBuild(kernels, nothing)
    @test backend_result.kernels === kernels
    @test backend_result.diagnostics === nothing
    @test QCLNEGF.validate_scattering_kernel_set(kernels, numerical, [acoustic]) ===
          kernels

    bad = KernelSet(
        Dict(:acoustic => normalized[:, 1:1, :, :, :, :]),
        Dict(:acoustic => 2.0),
        form_factor,
        [:acoustic],
    )
    @test_throws QCLNEGF.ScatteringValidationError begin
        QCLNEGF.validate_scattering_kernel_set(bad, numerical, [acoustic])
    end
end

end # independent suite
