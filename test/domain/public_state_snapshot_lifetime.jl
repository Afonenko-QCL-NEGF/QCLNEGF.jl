module PublicStateSnapshotLifetimeTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")

@testset "Exact snapshot declares its memory burst before copying accepted buffers" begin
    f=public_state_fixture()
    phases=NamedTuple[]
    snapshots=ExactSolverState[]
    production=QCLNEGF.with_production_options(f.production;phase_request=request->begin
        push!(phases,request)
        request.task_width
    end)
    N=QCLNEGF.QCLNumerics
    ctx=N._exact_context(f.problem,f.options,production,f.execution,f.policy,nothing,
        state->push!(snapshots,state))
    N._capture_inner(ctx,f.problem,f.state.published_U,f.state.scba,f.state.outer_history,
        Dict{String,Any}[],Dict{String,Any}[],1,0,f.state.continuation,0,0,f.state.mixer_position)
    @test length(phases)==1
    if !isempty(phases)
        @test phases[1].name===:exact_state_snapshot
        @test phases[1].stage===:checkpoint
        @test phases[1].memory_burst_bytes>=Base.summarysize(f.state.scba)
    end
    @test length(snapshots)==1
    if !isempty(snapshots)
        value=snapshots[1].scba.green.Gᴿ[1]
        f.state.scba.green.Gᴿ[1]=99.0+3im
        @test snapshots[1].scba.green.Gᴿ[1]==value
        pending=snapshots[1].continuation.candidate.plus.Σᴿ[1]
        f.state.continuation.candidate.plus.Σᴿ[1]=88.0+2im
        @test snapshots[1].continuation.candidate.plus.Σᴿ[1]==pending
    end
end
end
