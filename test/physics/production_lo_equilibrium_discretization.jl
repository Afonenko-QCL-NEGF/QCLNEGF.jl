module ProductionLOEquilibriumDiscretization
include("../support/common.jl")
const BN=QCLNEGF

function actual_lo_equilibrium(nodes, steps_per_phonon, algorithm)
    physical=reference_parameters(Tᴸ = 70u"K", V_period = 0u"mV")
    phonon=ustrip(u"eV", physical.ħωᴸᴼ)
    step=phonon/steps_per_phonon
    half_window=(nodes-1)*step/2
    original=tutorial_numerics()
    changes=(;
        N_b = 2,
        N_k = 3,
        N_E = nodes,
        E_min = -half_window*u"eV",
        E_max = half_window*u"eV",
    )
    numerical=NumericalParameters(;
        (
            name=>hasproperty(changes, name) ? getproperty(changes, name) :
                  getfield(original, name) for name in fieldnames(NumericalParameters)
        )...,
    )
    algorithms=AlgorithmOptions(energy_shift = algorithm)
    options=ProductionOptions(algorithms = algorithms, worker_count = 1)
    built=build_configured_scattering_problem(;
        physical,
        numerical,
        scales = ScaleSystem(),
        scattering = ScatteringOptions(
            LO = true,
            acoustic = false,
            impurity = false,
            IFR = false,
            alloy = false,
        ),
        algorithms,
        kernel_options = ProductionKernelOptions(),
    )
    problem=built.problem
    g=problem.grids
    s=problem.scales
    shape=(nodes, numerical.N_k, numerical.N_b, numerical.N_b)
    GR=zeros(ComplexF64, shape)
    less=similar(GR)
    greater=similar(GR)
    A=similar(GR)
    fill!(less, 0)
    fill!(greater, 0)
    fill!(A, 0)
    kBT=BN._electronvolts(CODATA.kᴮₑᵥ*physical.Tᴸ)/s.E₀_eV
    f=[BN._fermi(e, 0.0, kBT) for e in g.ε]
    # A positive smooth spectral fixture isolates the actual LO scattering
    # operator from a separate Dyson fixed point. The broad window suppresses
    # exterior-tail contamination, while the fractional shifts remain real.
    for e = 1:nodes, m = 1:numerical.N_k, a = 1:numerical.N_b
        spectral=exp(-((g.ε[e]*s.E₀_eV)/0.035)^2)*(1+0.1m+0.05a)
        A[e, m, a, a]=spectral
        GR[e, m, a, a]=-0.5im*spectral
        less[e, m, a, a]=im*f[e]*spectral
        greater[e, m, a, a]=-im*(1-f[e])*spectral
    end
    green=GreenState(
        GR,
        less,
        greater,
        A,
        ones(nodes, numerical.N_k),
        ones(nodes, numerical.N_k),
    )
    cache=build_production_cache(problem; options)
    lesser, greater=BN._production_lo_contraction(
        cache.kernels[:LO],
        green,
        problem,
        cache,
        problem.kernels.qᴷ[:LO],
        options,
    )
    kms=copy(lesser)
    number, energy, event_scale=0.0, 0.0, 0.0
    for e = 1:nodes, m = 1:numerical.N_k
        kms[e, m, :, :].-=f[e] .* (lesser[e, m, :, :] .- greater[e, m, :, :])
        incoming=real(
            tr(Matrix(view(lesser, e, m, :, :))*Matrix(view(green.Gᵍ, e, m, :, :))),
        )
        outgoing=real(
            tr(Matrix(view(greater, e, m, :, :))*Matrix(view(green.Gˡ, e, m, :, :))),
        )
        weight=g.wᴱ[e]*g.wᵏ[m]
        number+=weight*(incoming-outgoing)
        energy+=weight*g.ε[e]*(incoming-outgoing)
        event_scale+=weight*(abs(incoming)+abs(outgoing))
    end
    kms_residual=norm(kms)/max(norm(lesser), norm(greater))
    delta=phonon/s.E₀_eV
    return (;
        kms = kms_residual,
        number = abs(number)/event_scale,
        energy = abs(energy)/(delta*event_scale),
    )
end

@testset "Actual LO kernel: integer shifts, fractional KMS error and refinement" begin
    for algorithm in (:sparse_plan, :conservative_pair)
        integer=actual_lo_equilibrium(129, 8.0, algorithm)
        @test integer.kms<1e-12
        @test integer.number<1e-12
        @test integer.energy<1e-12
        coarse=actual_lo_equilibrium(129, 8.3, algorithm)
        fine=actual_lo_equilibrium(257, 16.6, algorithm)
        finer=actual_lo_equilibrium(513, 33.2, algorithm)
        # Finite-volume adjointness alone does not grant exact detailed balance
        # for an interpolated Fermi distribution at a fractional shift.
        @test coarse.kms>1e-8
        @test finer.kms<fine.kms<coarse.kms
        @test finer.energy<fine.energy<coarse.energy
        @test finer.number<1e-10
    end
end
end
