module FinalObservationContractTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")
const N=QCLNEGF.QCLNumerics

# Constructor-only failed terminal; no SCBA/Poisson or selected-map evaluation.
function final_fixture(;basis_states=1)
    f=public_state_fixture()
    problem=if basis_states==1
        f.problem
    else
        numerical=NumericalParameters(;(name=>(name===:N_b ? basis_states :
            getfield(f.problem.numerical,name)) for name in fieldnames(NumericalParameters))...)
        build_problem(;physical=f.problem.physical,numerical,
            scattering=ScatteringOptions(LO=false,acoustic=false,impurity=false,IFR=false,alloy=false))
    end
    ne,nk,nb=problem.numerical.N_E,problem.numerical.N_k,problem.numerical.N_b
    shape=(ne,nk,nb,nb)
    green=GreenState(fill(1.0-1im,shape),fill(0.5im,shape),fill(-1.5im,shape),
        fill(2.0+0im,shape),ones(ne,nk),ones(ne,nk))
    family()=SelfEnergyFamily(zeros(ComplexF64,shape),zeros(ComplexF64,shape),zeros(ComplexF64,shape))
    inner=public_state_failed_scba(f.state.scba)
    scba=SCBAResult(green,Dict{Symbol,SelfEnergyFamily}(),family(),family(),family(),
        inner.history,false,:invalid_candidate,:invalid,inner.restart_contract,inner.mixer_state)
    observations=Dict{Symbol,Any}(:quality=>:invalid,:warnings=>
        [N._core_plain(row) for row in f.state.warnings],:acceptance_measurement_context=>Dict(
            "backend"=>"production","hilbert_roundoff_measured"=>true,
            "map_kernels_shared_with_solver"=>true,"final_state_alignment"=>"producer_aligned"))
    report=ConvergenceReport(false,Dict(:r_K=>NaN,:r_U=>Inf),["constructor-only failed evidence"])
    solution=NEGFSolution(problem,f.options,zeros(problem.numerical.N_z),ones(problem.numerical.N_z),
        scba,OuterIteration[],observations,report,false,:scba_invalid_candidate)
    identity=exact_restart_identity(problem,f.options,f.production,f.execution;
        domain_adaptation=f.policy,kernel_options=nothing)
    state=replace_public_state(f.state;identity,boundary=:final,scba,
        published_U=solution.Uᴴ,published_n=solution.n,continuation=nothing,
        status=solution.status,final_evidence=N._final_evidence(solution))
    binding=(;options=f.options,production_options=f.production,execution=f.execution,
        domain_adaptation=f.policy,kernel_options=nothing)
    validate_exact_state(problem,state;binding...)
    return (;f...,problem,solution,state,binding)
end
function replace_solution(s;changes...)
    names=fieldnames(NEGFSolution)
    args=merge(NamedTuple{names}(Tuple(getfield(s,n) for n in names)),(;changes...))
    NEGFSolution((getfield(args,n) for n in names)...)
end
function take_snapshot(f;source_mode=:produced,kwargs...)
    final_observation_snapshot(f.solution,f.state;f.binding...,source_mode,
        maximum_bytes=10^8,maximum_metadata_bytes=10^7,kwargs...)
end
caught(call)=try call();nothing catch error;error end

@testset "Final observation public port is defined" begin
    for symbol in (:FinalObservationSnapshot,:FinalObservationMemoryEstimate,
        :final_observation_memory_estimate,:final_observation_snapshot)
        @test isdefined(QCLNEGF,symbol)
    end
end

if isdefined(QCLNEGF,:final_observation_snapshot)
@testset "Same-call terminal binding refuses swaps before observations" begin
    f=final_fixture()
    @test take_snapshot(f).terminal_integrity==f.state.integrity
    for key in (:Uᴴ,:n)
        array=copy(getfield(f.solution,key));array[1]+=1
        wrong=replace_solution(f.solution;NamedTuple{(key,)}((array,))...)
        err=caught(()->take_snapshot((;f...,solution=wrong)))
        @test err isa ArgumentError
        @test occursin("binding",sprint(showerror,err))
    end
    for name in fieldnames(GreenState)
        wrong=deepcopy(f.solution)
        getfield(wrong.scba.green,name)[1]+=1
        @test caught(()->take_snapshot((;f...,solution=wrong))) isa ArgumentError
    end
    for group in (:embedding,:embedding_plus,:embedding_minus), name in (:Σᴿ,:Σˡ,:Σᵍ)
        wrong=deepcopy(f.solution);getfield(getfield(wrong.scba,group),name)[1]+=1
        @test caught(()->take_snapshot((;f...,solution=wrong))) isa ArgumentError
    end
    for mutate in (s->push!(s.scba.history,s.scba.history[1]),
        s->push!(s.scba.mixer_state.states,[s.scba.embedding_plus]),
        s->(s.report.metrics[:r_K]=0.0),
        s->(s.observables[:warnings][1]["code"]="foreign warning"),
        s->(s.observables[:quality]=:strict))
        wrong=deepcopy(f.solution);mutate(wrong)
        @test caught(()->take_snapshot((;f...,solution=wrong))) isa ArgumentError
    end
    options=SolverOptions(;(name=>(name===:max_scba ? 4 : getfield(f.options,name))
        for name in fieldnames(SolverOptions))...)
    @test caught(()->take_snapshot((;f...,solution=replace_solution(f.solution;options)))) isa ArgumentError
    @test caught(()->take_snapshot(f;source_mode=:guessed)) isa ArgumentError
    estimate=final_observation_memory_estimate(f.solution,f.state;maximum_metadata_bytes=10^7)
    @test estimate.peak_additional_bytes>estimate.view_bytes>0
    @test estimate.caller_resident_bytes>0
    @test caught(()->final_observation_snapshot(f.solution,f.state;f.binding...,
        source_mode=:produced,maximum_bytes=1,maximum_metadata_bytes=10^7)) isa ArgumentError
    @test caught(()->final_observation_memory_estimate(f.solution,f.state;
        maximum_metadata_bytes=1)) isa ArgumentError
end

function replace_observation(value;changes...)
    names=fieldnames(typeof(value))
    args=merge(NamedTuple{names}(Tuple(getfield(value,n) for n in names)),(;changes...))
    typeof(value)((getfield(args,n) for n in names)...)
end
@testset "Public observation constructors enforce structural links, not new physics gates" begin
    f=final_fixture();snap=take_snapshot(f);c=snap.context
    @test_throws DimensionMismatch replace_observation(c;wE=Float64[])
    @test_throws ArgumentError replace_observation(c;domain_revision=-1)
    @test_throws ArgumentError replace_observation(c;state_index=Int64[99])
    # Serialization can retain failed/raw weights and limiting temperatures.
    @test replace_observation(c;wE=fill(-1.0,length(c.wE))).wE[1]==-1.0
    inputs=replace_observation(c.inputs;T_LO_K=0.0)
    @test replace_observation(c;inputs).inputs.T_LO_K==0.0
    @test_throws DimensionMismatch replace_observation(snap.views;density_per_m3=Float64[])
    @test_throws ArgumentError replace_observation(snap;source_mode=:invented)
    @test_throws ArgumentError replace_observation(snap;positivity=nothing)
    @test_throws OverflowError N._final_numeric_sizes((typemax(Int),2,2,2))
    @test N._final_numeric_sizes((5,3,2,2)).context>=264 # one UInt64 BitVector word, not three bytes
end
@testset "LAPACK workspace forecast queries scalar buffers, not an eigensystem" begin
    @test isdefined(N,:_final_eigen_workspace_bytes)
    if isdefined(N,:_final_eigen_workspace_bytes)
        one=N._final_eigen_workspace_bytes(1)
        large=N._final_eigen_workspace_bytes(100_000)
        @test large>one>0
        # 100000² complex values would require 160 GB; the metadata query owns scalar refs only.
        allocated=@allocated N._final_eigen_workspace_bytes(100_000)
        @test allocated<65536
        @test_throws ArgumentError N._final_eigen_workspace_bytes(0)
        @test_throws OverflowError N._final_eigen_workspace_bytes(typemax(Int))
        @test_throws ArgumentError N._final_query_count(Inf)
        @test_throws ArgumentError N._final_query_count(1.5)
        f=final_fixture()
        base=N._final_numeric_sizes((25,9,2,1)).workspace
        actual=final_observation_memory_estimate(f.solution,f.state;maximum_metadata_bytes=10^7)
        @test actual.workspace_bytes==base+one
    end
end
end
end
