module PublicStateContractTests
println("[P02] loading prepared Core for public state constructors")
flush(stdout)
include("../support/common.jl")

@testset "Public exact-state constructors reject ambiguous identity" begin
    @test isdefined(QCLNEGF, :ExactExecutionIdentity)
    @test isdefined(QCLNEGF, :exact_restart_identity)
    if isdefined(QCLNEGF, :ExactExecutionIdentity)
        @test_throws ArgumentError ExactExecutionIdentity(
            core_code="", julia="1.13.0", blas="fixture-blas", fftw="fixture-fftw",
            worker_count=1, blas_threads=1, energy_chunk=8, hilbert_columns=2,
            residual_chunk=8, parallel_backend=:threads, scalar_type="Float64/ComplexF64",
            index_type="Int64", reduction="static_chunks_v1",
        )
    end
end

@testset "Scalar events are immutable values with explicit absence" begin
    @test isdefined(QCLNEGF, :ScientificHistorySnapshot)
    metric=SolverMetric(:current, 0.0, "A/m^2")
    event=SolverEvent(:update, :scba, "accepted", :running, 1, 2, [metric], "")
    @test event.metrics isa Tuple
    @test metric.value === 0.0
    if hasproperty(metric, :measurement_status)
        @test metric.measurement_status === :measured
        absent=SolverMetric(:previous_current, nothing, "A/m^2";
            measurement_status=:not_measured, reason=:no_previous_iteration)
        @test absent.value === nothing
        @test absent.reason === :no_previous_iteration
        @test_throws ArgumentError SolverMetric(:current, nothing, "A/m^2";
            measurement_status=:measured)
        nonfinite=SolverMetric(:raw_residual,NaN,"1")
        @test isnan(nonfinite.value)
        @test nonfinite.measurement_status===:error
        @test nonfinite.reason===:non_finite
    else
        @test false
    end
end

@testset "Scientific row transfer detaches every nested array" begin
    if isdefined(QCLNEGF, :ScientificHistorySnapshot)
        populations=[0.25,0.75]
        density=[1.0,2.0,3.0]
        row=OuterIteration(1, 0.1,0.2,0.3,0.0,0.0,Inf,Inf,0.0,0.0,
            populations,density)
        snapshot=ScientificHistorySnapshot(stage=:outer, outer_iteration=1,
            domain_revision=0,discretization_id="explicit grid",basis_id="explicit basis",
            inner_iteration=nothing, sequence=1, converged=false,
            status=:running, quality=:unresolved, strict_streak=0,
            diagnostic_streak=0, row=row)
        populations[1]=9.0
        density[1]=9.0
        @test snapshot.row.populations == [0.25,0.75]
        @test snapshot.row.density == [1.0,2.0,3.0]
        snapshot.row.density[2]=8.0
        @test row.density[2] == 2.0
        @test snapshot.domain_revision==0 && snapshot.discretization_id=="explicit grid" &&
            snapshot.basis_id=="explicit basis"
        @test hasproperty(snapshot,:missing_reasons)
        @test hasproperty(snapshot,:units)
        if hasproperty(snapshot,:missing_reasons) && hasproperty(snapshot,:units)
            @test Dict(snapshot.missing_reasons)[:r_Jchange]===:no_previous_iteration
            @test !haskey(Dict(snapshot.missing_reasons),:J)
            @test Dict(snapshot.units)[:J]=="A/m^2"
            @test Dict(snapshot.units)[:density]=="L0^3 n"
        end
    else
        @test false
    end
end

@testset "Canonical identity uses fixed binary bits rather than approximate values" begin
    if isdefined(QCLNEGF.QCLNumerics, :_exact_canonical_bytes)
        bytes=QCLNEGF.QCLNumerics._exact_canonical_bytes
        # v1 scalar tag f + IEEE754 binary64, network byte order.
        @test bytes(1.0) == UInt8[0x66,0x3f,0xf0,0,0,0,0,0,0]
        @test bytes(-0.0) == UInt8[0x66,0x80,0,0,0,0,0,0,0]
        @test bytes(0.0) != bytes(-0.0)
        @test bytes(1.0) != bytes(nextfloat(1.0))
        @test bytes2hex(UInt8[QCLNEGF.QCLNumerics._exact_digest(1.0)...]) ==
            "e2ac13292d1245bac0ce7b682e9b0efef04e53d29a40be9d72eb2b6bc21b3533"
        @test bytes(Dict("b"=>2,"a"=>1)) == bytes(Dict("a"=>1,"b"=>2))
        @test bytes(reshape([1.0,2.0],2,1)) != bytes(reshape([1.0,2.0],1,2))
    else
        @test false
    end
end

@testset "Scoped evidence snapshots reject ambiguous keys and isolate arrays" begin
    if isdefined(QCLNEGF, :CoreRecord)
        source=Dict("metric"=>[1.0,2.0],"status"=>"not_measured")
        record=CoreRecord(source)
        source["metric"][1]=9.0
        @test record["metric"].values == [1.0,2.0]
        @test record["status"].value == "not_measured"
        @test_throws ArgumentError CoreRecord(("a"=>CoreScalar(1),"a"=>CoreScalar(2)))
        @test_throws ArgumentError CoreRecord(Dict("bad"=>sin))
    else
        @test false
    end
end

@testset "Domain revision variants declare their required source fields" begin
    @test isdefined(QCLNEGF,:DomainPreflightEvidence)
    @test isdefined(QCLNEGF,:DomainMeasuredEvidence)
    if isdefined(QCLNEGF,:DomainPreflightEvidence) && isdefined(QCLNEGF,:DomainMeasuredEvidence)
        before=DomainPreflightEvidence(discretization_id="d",basis_id="b",E_min_eV=-1.0,
            E_max_eV=1.0,NE=3,energy_min_scaled=-1.0,energy_max_scaled=1.0,
            required_min_scaled=-2.0,required_max_scaled=2.0,
            reason="estimate",seed_transfer="none; no state solved")
        row=DomainRevision(revision=0,source=:preflight_estimate,evidence=before)
        @test row.evidence isa DomainPreflightEvidence
        @test_throws ArgumentError DomainRevision(revision=0,source=:measured_state,evidence=before)
        @test_throws UndefKeywordError DomainPreflightEvidence(discretization_id="d")
        metrics=Dict(:r_K=>NaN)
        measured=DomainMeasuredEvidence(discretization_id="d",basis_id="b",hartree_id="h",
            E_min_eV=-1.0,E_max_eV=1.0,NE=3,energy_min_scaled=-1.0,energy_max_scaled=1.0,
            required_min_scaled=-2.0,required_max_scaled=2.0,tail_low=NaN,tail_high=Inf,
            edge_gamma=-Inf,edge_spectral=0.0,scba_quality=:invalid,
            seconds=1.0,status="failed",strict_accepted=false,
            metrics=metrics,current_A_per_m2=Inf,momentum_tail_unresolved=true,
            seed_transfer="cold rebuild",reason="measured")
        metrics[:r_K]=0.0
        @test isnan(measured.metrics[:r_K]) && isinf(measured.current_A_per_m2)
        @test isnan(measured.tail_low) && measured.tail_high===Inf && measured.edge_gamma===-Inf
        @test_throws UndefKeywordError DomainMeasuredEvidence(discretization_id="d")
        for (source,evidence) in ((:preflight_estimate,before),(:measured_state,measured))
            complete=QCLNEGF.QCLNumerics._domain_revision_plain(
                DomainRevision(;revision=0,source,evidence))
            for name in fieldnames(typeof(evidence))
                missing=copy(complete)
                delete!(missing,String(name))
                @test_throws ArgumentError QCLNEGF.QCLNumerics._domain_revisions([missing])
            end
        end
        legacy=QCLNEGF.QCLNumerics._domain_revision_plain(
            DomainRevision(revision=0,source=:measured_state,evidence=measured))
        @test isequal(only(QCLNEGF.QCLNumerics._domain_revisions([legacy])).evidence.metrics,measured.metrics)
        delete!(legacy,"hartree_id")
        @test_throws ArgumentError QCLNEGF.QCLNumerics._domain_revisions([legacy])
    end
end
end
