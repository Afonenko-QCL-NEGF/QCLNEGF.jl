module Suite_T007
include("../support/common.jl")
include("../support/algorithm_modes.jl")

@testset "Anderson dependent histories and safeguards" begin
    family(value) = SelfEnergyFamily(
        fill(ComplexF64(value), 1, 1, 1, 1),
        zeros(ComplexF64, 1, 1, 1, 1),
        zeros(ComplexF64, 1, 1, 1, 1),
    )
    algorithms = AlgorithmOptions(
        mixing = :anderson,
        anderson_history_depth = 2,
        anderson_regularization = 0.0,
    )
    options = ProductionOptions(algorithms = algorithms)
    workspace = QCLNEGF._AndersonWorkspace()

    QCLNEGF._anderson_mix!(workspace, [family(0)], [family(1)], algorithms, 1.0, options)
    repeated = [family(0)]
    QCLNEGF._anderson_mix!(workspace, repeated, [family(1)], algorithms, 1.0, options)
    # Duplicate residuals are a valid rank-one history. The declared
    # rank-revealing solve must produce the same bounded fixed-point step.
    coefficients = QCLNEGF._anderson_coefficients(workspace.residuals, 0.0)
    @test coefficients ≈ [0.5, 0.5]
    for component in (:Σᴿ, :Σˡ, :Σᵍ)
        @test getfield(only(repeated), component) ≈ getfield(family(1.0), component)
    end

    # A fivefold residual increase exceeds the squared-norm safeguard and
    # discards the old history before taking the damped current step.
    growing = [family(0)]
    QCLNEGF._anderson_mix!(workspace, growing, [family(5)], algorithms, 0.25, options)
    @test length(workspace.states) == length(workspace.residuals) == 1
    for component in (:Σᴿ, :Σˡ, :Σᵍ)
        @test getfield(only(growing), component) ≈ getfield(family(1.25), component)
    end

    # At an exact fixed point, the zero Gram system has no usable history;
    # the specified linear fallback must preserve the fixed point exactly.
    fixed_workspace = QCLNEGF._AndersonWorkspace()
    for _ = 1:2
        fixed = [family(2)]
        QCLNEGF._anderson_mix!(
            fixed_workspace,
            fixed,
            [family(2)],
            algorithms,
            1.0,
            options,
        )
        for component in (:Σᴿ, :Σˡ, :Σᵍ)
            @test getfield(only(fixed), component) == getfield(family(2.0), component)
        end
    end
    @test_throws ArgumentError QCLNEGF._anderson_coefficients(
        fixed_workspace.residuals,
        0.0,
    )

    # Even a modest resolved deterioration resets extrapolation history.
    modest = QCLNEGF._AndersonWorkspace()
    QCLNEGF._anderson_mix!(modest, [family(0)], [family(1)], algorithms, 0.25, options)
    output=[family(0)]
    QCLNEGF._anderson_mix!(modest, output, [family(1.1)], algorithms, 0.25, options)
    @test length(modest.states)==length(modest.residuals)==1
    @test only(output).Σᴿ ≈ fill(0.275+0im, 1, 1, 1, 1)
end

end # independent suite
