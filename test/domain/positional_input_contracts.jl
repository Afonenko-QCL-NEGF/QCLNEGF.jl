module PositionalInputContracts
include("../support/common.jl")
function positional(T, value, field, replacement)
    T((name===field ? replacement : getfield(value, name) for name in fieldnames(T))...)
end
@testset "Positional public physical and numerical constructors enforce inputs" begin
    p=reference_parameters()
    n=tutorial_numerics()
    scales=ScaleSystem()
    opts=SolverOptions()
    @test_throws ArgumentError positional(PhysicalParameters, p, :Tᴸ, -1u"K")
    @test_throws ArgumentError positional(PhysicalParameters, p, :f_ion, 2.0)
    @test_throws ArgumentError positional(PhysicalParameters, p, :F_bias, NaN*u"V/m")
    @test_throws ArgumentError positional(NumericalParameters, n, :N_E, 1)
    @test_throws ArgumentError positional(NumericalParameters, n, :E_max, n.E_min)
    @test_throws ArgumentError positional(ScaleSystem, scales, :E₀_eV, 2scales.E₀_eV)
    @test_throws ArgumentError positional(SolverOptions, opts, :α_Σ, 0.0)
    @test_throws ArgumentError positional(SolverOptions, opts, :max_scba, 0)
    @test_throws ArgumentError Layer(-1u"nm", :GaAs, 0.0, 0u"eV", 0.067, 0.067, 12.9, false)
end
end
