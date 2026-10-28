module PublicOptionsContractTests
include("../support/common.jl")
const N=QCLNEGF.QCLNumerics

@testset "Public option validation rejects invalid families without executing work" begin
    available=isdefined(QCLNEGF,:validate_options) && isdefined(QCLNEGF,:validate_result_classification)
    @test available
    if available
        called=Ref(false)
        production=ProductionOptions(worker_count=1,event_sink=_->(called[]=true),
            phase_request=_->(called[]=true))
        for value in (SolverOptions(),AlgorithmOptions(),production,ProductionKernelOptions())
            @test validate_options(value)===value
        end
        for make in (() -> SolverOptions(α_Σ=0.0),() -> AlgorithmOptions(energy_shift=:cyclic),
            () -> ProductionOptions(worker_count=-1),() -> ProductionKernelOptions(relative_tolerance=1.0))
            @test_throws ArgumentError validate_options(make())
        end
        @test !called[]
        @test validate_result_classification(true,:converged,:strictly_converged)===:strictly_converged
        @test validate_result_classification(false,:approximate,:strictly_converged)===:strictly_converged
        @test validate_result_classification(false,:failed,:invalid)===:invalid
        @test_throws ArgumentError validate_result_classification(false,:converged,:strictly_converged)
        @test_throws ArgumentError validate_result_classification(true,:converged,:unresolved)
        @test_throws ArgumentError validate_result_classification(false,:failed,:unknown)
    end
end

@testset "Prebuild memory classification is public but its estimate does not grant a budget" begin
    n=tutorial_numerics()
    available=hasmethod(estimate_production_memory,Tuple{NumericalParameters,ScatteringOptions})
    @test available
    if available
        production=ProductionOptions(worker_count=1,memory_budget_bytes=1)
        # Literal mechanism/dense/LO classifications, independent of the new classifier.
        cases=((ScatteringOptions(LO=false,acoustic=false,impurity=false,IFR=false,alloy=false),0,0,:absent),
            (ScatteringOptions(LO=false,acoustic=true,impurity=false,IFR=false,alloy=true),2,0,:absent),
            (ScatteringOptions(LO=true,acoustic=false,impurity=false,IFR=false,alloy=false),1,1,:dense),
            (ScatteringOptions(LO=true,acoustic=true,impurity=true,IFR=true,alloy=true),5,3,:dense))
        for (scattering,count,dense,lo) in cases
            estimate=estimate_production_memory(n,scattering;options=production,worker_capacity=1)
            literal=estimate_production_memory(n,count;dense_mechanism_count=dense,lo_kernel=lo,
                options=production,worker_capacity=1)
            @test estimate.peak_bytes==literal.peak_bytes>production.memory_budget_bytes
            @test estimate.breakdown==literal.breakdown
            # Existing admission remains enforced; do not bypass it by changing options.
            @test_throws ArgumentError N._preflight_production_memory(n,scattering,production)
            ample=ProductionOptions(worker_count=1,memory_budget_bytes=typemax(Int))
            @test N._preflight_production_memory(n,scattering,ample).peak_bytes==
                estimate_production_memory(n,scattering;options=ample).peak_bytes
        end
    end
end
end
