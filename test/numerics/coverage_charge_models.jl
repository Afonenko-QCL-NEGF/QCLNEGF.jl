module CoverageChargeModelContracts
include("../support/common.jl")

@testset "Independent coverage, particle constraint and optional physical models" begin
    charge = charge_constraint_diagnostics(1.2, 1.0, 1.0; previous_lambda = 1.3)
    @test charge.constraint_residual ≈ 0.2
    @test charge.normalization_residual == 0
    @test charge.lambda_change ≈ -0.1
    @test_throws ArgumentError charge_constraint_diagnostics(NaN, 1, 1)
    @test_throws ArgumentError PhysicalModelOptions(lo_population = :rate_balance)
    @test_throws ArgumentError PhysicalModelOptions(screening = :density_temperature)
    physical = reference_parameters(Tᴸ = 70u"K")
    no_scattering = ScatteringOptions(
        LO = false,
        acoustic = false,
        impurity = false,
        IFR = false,
        alloy = false,
    )
    problem = build_problem(;
        physical,
        numerical = tutorial_numerics(),
        scattering = no_scattering,
    )
    coverage = representation_coverage(problem)
    h = project_hamiltonians(problem, zeros(problem.numerical.N_z))
    direct = [eigvals(Hermitian(Matrix(view(h, m, :, :)))) for m in axes(h, 1)]
    @test coverage.h_min ≈ minimum(minimum.(direct))
    @test coverage.h_max ≈ maximum(maximum.(direct))
    @test coverage.source === :hamiltonian_estimate
    @test coverage.discretization_id == representation_coverage(problem).discretization_id
    shifted = representation_coverage(problem; hartree = fill(0.2, problem.numerical.N_z))
    @test shifted.h_max ≈ coverage.h_max+0.2 atol=1e-12
    @test shifted.hartree_id != coverage.hartree_id
    @test shifted.basis_id == coverage.basis_id
    expanded = expand_energy_window(problem, coverage)
    @test expanded.E_max >= problem.numerical.E_max
    @test (expanded.E_max-expanded.E_min)/(expanded.N_E-1) <=
          (problem.numerical.E_max-problem.numerical.E_min)/(problem.numerical.N_E-1)*(
        1+1e-14
    )
    # Lorentzian oracle integrates independently by the trapezoid rule.
    e = collect(range(-2.0, 3.0; length = 20_001))
    step = e[2]-e[1]
    spectrum = @. 0.1/(π*((e-0.3)^2+0.1^2))
    quadrature = step*(sum(spectrum)-(first(spectrum)+last(spectrum))/2)
    @test quadrature ≈ lorentzian_window_weight(-2.0, 3.0, 0.3, 0.1) atol=1e-10
    @test lorentzian_window_weight(-20.0, 30.0, 0.3, 0.1) > quadrature
    # A coherent local transverse dispersion modifies Hamiltonian, not dk weights.
    kane = PhysicalModelOptions(dispersion = :kane_inplane, nonparabolicity_per_eV = 1.2)
    variant = NEGFProblem(
        problem.physical,
        problem.numerical,
        problem.scattering,
        problem.scales,
        problem.grids,
        problem.profiles,
        problem.basis,
        problem.kernels,
        problem.W₊ᴱᵖ,
        problem.W₋ᴱᵖ,
        problem.W₊ᴸᴼ,
        problem.W₋ᴸᴼ,
        problem.energy_shift_discretization,
        kane,
    )
    hn = project_hamiltonians(variant, zeros(problem.numerical.N_z))
    @test hn[1, :, :] ≈ h[1, :, :]
    @test eigmax(Hermitian(Matrix(hn[end, :, :]-h[end, :, :]))) <= 1e-12
    @test variant.grids.wᵏ == problem.grids.wᵏ
    for energy in (0.0, 0.1, 3.0)
        transformed = QCLNEGF.transverse_kinetic_energy(energy, 0.8)
        @test transformed*(1+0.8transformed) ≈ energy atol=1e-15
    end
    screened = QCLNEGF.resolve_physical_models(
        physical,
        PhysicalModelOptions(
            screening = :density_temperature,
            screening_temperature_K = 70,
            screening_mass_ratio = 0.067,
        ),
    )
    @test screened.q_s > 0u"m^-1" && screened.qᴸᴼ_s > 0u"m^-1"
    @test screened.layers == physical.layers
    @test screened.q_s != physical.q_s
end
end
