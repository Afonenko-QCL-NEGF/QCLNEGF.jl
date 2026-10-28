module PublicStateRejectionTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")

@testset "Incomplete or incompatible exact groups fail before numerical execution" begin
    f=public_state_fixture(scattering_options=ScatteringOptions(
        LO=false,acoustic=true,impurity=false,IFR=false,alloy=false))
    bind=(options=f.options,production_options=f.production,execution=f.execution,
        domain_adaptation=f.policy,kernel_options=nothing)
    @test validate_exact_state(f.problem,f.state;bind...)===nothing
    state=deepcopy(f.state)
    delete!(state.scba.scattering,:acoustic)
    @test_throws ArgumentError validate_exact_state(f.problem,replace_public_state(state);bind...)
    for field in (:Σᴿ,:Σˡ,:Σᵍ), group in (:scattering,:embedding,:embedding_plus,:embedding_minus,:candidate)
        state=deepcopy(f.state)
        family=group===:scattering ? state.scba.scattering[:acoustic] :
            group===:candidate ? state.continuation.candidate.scattering[:acoustic] : getfield(state.scba,group)
        getfield(family,field)[1]=ComplexF64(NaN,0)
        @test_throws ArgumentError validate_exact_state(f.problem,replace_public_state(state);bind...)
    end
    wrong=replace_public_state(f.state;family_order=(:embedding_plus,:acoustic,:embedding_minus))
    @test_throws ArgumentError validate_exact_state(f.problem,wrong;bind...)
    wrong=replace_public_state(f.state;progress=SolverProgress(0,1,2,0,0,0))
    @test_throws ArgumentError validate_exact_state(f.problem,wrong;bind...)
    wrong=replace_public_state(f.state;continuation=nothing)
    @test_throws ArgumentError validate_exact_state(f.problem,wrong;bind...)
    wrong=replace_public_state(f.state;boundary=:final,continuation=nothing,final_evidence=nothing)
    @test_throws ArgumentError validate_exact_state(f.problem,wrong;bind...)
    v=f.state.identity
    wrong=replace_public_state(f.state;identity=ExactRestartIdentity(UInt16(2),v.problem,v.policy,v.execution))
    @test_throws ArgumentError validate_exact_state(f.problem,wrong;bind...)
    old=f.state.scba
    families=[deepcopy(old.scattering[:acoustic]),deepcopy(old.embedding_plus),deepcopy(old.embedding_minus)]
    mixer=SCBAMixerState(:linear,[families],[deepcopy(families)])
    scba=SCBAResult(old.green,old.scattering,old.embedding,old.embedding_plus,old.embedding_minus,
        old.history,old.converged,old.status,old.quality,old.restart_contract,mixer)
    wrong=replace_public_state(f.state;scba=scba)
    @test_throws ArgumentError validate_exact_state(f.problem,wrong;bind...)
    workers=Base.Threads.nthreads(:default)+1
    oversized=QCLNEGF.with_production_options(f.production;worker_count=workers)
    claimed=ExactExecutionIdentity(;
        (name=>(name===:worker_count ? workers : getfield(f.execution,name))
         for name in fieldnames(ExactExecutionIdentity))...)
    @test QCLNEGF.QCLNumerics._production_thread_worker_count(oversized,typemax(Int))==workers-1
    @test_throws ArgumentError exact_restart_identity(f.problem,f.options,oversized,claimed;
        domain_adaptation=f.policy,kernel_options=nothing)
end
end
