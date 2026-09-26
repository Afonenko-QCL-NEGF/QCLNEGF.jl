module PackageBoundaryTests
include("../support/common.jl")
@testset "Independent in-memory package" begin
    @test !isdefined(QCLNEGF, :YAML)
    @test !isdefined(QCLNEGF, :HDF5)
    @test !isdefined(QCLNEGF, :QCLNEGFRunner)
    @test !isdefined(QCLNEGF, :load_run_configuration)
    @test !isdefined(QCLNEGF, :save_checkpoint)
    @test parentmodule(QCLNEGF.solve_production) === QCLNEGF
    @test QCLNEGF._software_version() == "0.2.0"
end
end
