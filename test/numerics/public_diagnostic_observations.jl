module PublicDiagnosticObservationsTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")
const N=QCLNEGF.QCLNumerics

function manufactured_solution(;Nb=1)
    f=public_state_fixture()
    n=NumericalParameters(; (name=>(name===:N_b ? Nb : getfield(f.problem.numerical,name))
        for name in fieldnames(NumericalParameters))...)
    # No SCBA/Poisson: this is an explicit matrix-observation fixture.
    problem=build_problem(physical=f.problem.physical,numerical=n,scales=f.problem.scales,
        scattering=f.problem.scattering,physical_models=f.problem.models,validate_static=false)
    shape=(n.N_E,n.N_k,Nb,Nb)
    GR=zeros(ComplexF64,shape); GL=similar(GR); GG=similar(GR); A=similar(GR)
    fill!(GL,0);fill!(GG,0);fill!(A,0)
    for e in 1:n.N_E,k in 1:n.N_k
        if Nb==2
            GR[e,k,:,:]=ComplexF64[1-im 1;0 2-im]
            GL[e,k,:,:]=im*ComplexF64[2 im;-im 1]
            GG[e,k,:,:]=-im*Matrix{ComplexF64}(I,2,2)
            A[e,k,:,:]=ComplexF64[-2 0;0 1]
        else
            GR[e,k,1,1]=1-im;GL[e,k,1,1]=im;GG[e,k,1,1]=-im;A[e,k,1,1]=1
        end
    end
    family()=SelfEnergyFamily(zeros(ComplexF64,shape),zeros(ComplexF64,shape),zeros(ComplexF64,shape))
    embedding=family()
    for e in 1:n.N_E,k in 1:n.N_k
        if Nb==2
            embedding.Σˡ[e,k,:,:]=im*ComplexF64[-1 0;0 2]
            embedding.Σᵍ[e,k,:,:]=-im*Matrix{ComplexF64}(I,2,2)
        end
    end
    green=GreenState(GR,GL,GG,A,ones(n.N_E,n.N_k),ones(n.N_E,n.N_k))
    scba=SCBAResult(green,Dict{Symbol,SelfEnergyFamily}(),embedding,family(),family(),SCBAIteration[],
        false,:running,:unresolved,nothing,SCBAMixerState())
    metrics=Dict(:r_K=>0.25,:r_sum=>0.75,:r_C_LO=>NaN,:r_D_candidate=>0.5)
    solution=NEGFSolution(problem,f.options,zeros(n.N_z),zeros(n.N_z),scba,OuterIteration[],
        Dict{Symbol,Any}(),ConvergenceReport(false,metrics,["manufactured metadata"]),false,:unconverged)
    return solution
end

@testset "Metric diagnostics preserve final values, dynamic limits and explicit disabled approximation" begin
    available=isdefined(QCLNEGF,:solution_metric_diagnostics)
    @test available
    if available
        solution=manufactured_solution()
        result=solution_metric_diagnostics(solution)
        @test result.summary_metrics[:keldysh]==0.25 && result.summary_metrics[:spectral_sum]==0.75
        @test result.strict_limits[:r_C_LO]==solution.options.tolerances.r_C
        @test result.strict_limits[:r_D_candidate]==solution.options.tolerances.r_K
        @test result.approximate_limits===nothing
        @test isnan(solution.report.metrics[:r_C_LO])
        result.summary_metrics[:keldysh]=100.0
        @test solution.report.metrics[:r_K]==0.25
        @test !solution_scientific_assessment(solution)["scientific_accepted"]
        push!(solution.scba.history,SCBAIteration(1,0.2,0.3,0.4,0.5,0.0,1.0,0.0,0.0,0.0,
            Inf,Inf,3.0,1.0,1.0,1.0,NaN,nothing,nothing))
        push!(solution.outer_history,OuterIteration(1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
            [1.0],ones(solution.problem.numerical.N_z)))
        mixed=solution_metric_diagnostics(solution)
        @test mixed.summary_metrics[:keldysh]==0.25 # final report overrides history's 0.4
        @test mixed.summary_metrics[:poisson]==0.2 # history retained if report has no entry
        solution.report.metrics[:r_roundoff]=NaN
        @test isnan(solution_metric_diagnostics(solution).summary_metrics[:fft_roundoff])
        q=DiagnosticQualityPolicy(enabled=true)
        options=SolverOptions(convergence=ConvergencePolicy(diagnostic_quality=q))
        approximate=NEGFSolution((name===:options ? options : getfield(solution,name)
            for name in fieldnames(NEGFSolution))...)
        bands=solution_metric_diagnostics(approximate)
        @test bands.approximate_limits[:r_K]==q.keldysh_threshold
        @test bands.approximate_limits[:r_PSD]==q.positivity_threshold
        @test bands.strict_limits[:r_K]==options.tolerances.r_K
        @test bands.approximate_limits[:r_D]==bands.strict_limits[:r_D]
        @test bands.approximate_limits[:r_C_LO]==bands.strict_limits[:r_C_LO]
        bands.approximate_limits[:r_C_LO]=100
        @test solution_metric_diagnostics(approximate).approximate_limits[:r_C_LO]==options.tolerances.r_C
        @test !solution_scientific_assessment(approximate)["scientific_accepted"]
    end
end

@testset "Six matrix kinds retain full integrals and literal noncommuting raw Keldysh defects" begin
    available=isdefined(QCLNEGF,:positivity_diagnostics)
    @test available
    if available
        solution=manufactured_solution(Nb=2)
        problem=solution.problem; ne,nk=problem.numerical.N_E,problem.numerical.N_k
        light=positivity_diagnostics(solution)
        @test light.integrals===nothing && light.scan_order===:energy_momentum
        @test light.matrix_kinds==(:spectral,:occupied,:unoccupied,:broadening,:raw_occupied,:raw_unoccupied)
        full=positivity_diagnostics(solution;integrate=true)
        @test full.scan_order===:momentum_energy
        @test full.witnesses[1].minimum_eigenvalue≈-2.0
        @test full.witnesses[5].minimum_eigenvalue≈5-3sqrt(5.0) rtol=3e-15
        @test all(w->w.status===:measured,full.witnesses)
        @test full.integrals.negative_block_count==(ne*nk,0,0,0,ne*nk,0)
        @test full.integrals.invalid_block_count==(0,0,0,0,0,0)
        @test !full.strict_threshold_satisfied[1] && !full.strict_threshold_satisfied[5]
        energy_weight=sum(problem.grids.wᴱ)/(2π)
        for k in 1:nk
            @test full.integrals.spectral_integral[k,:,:]≈energy_weight*ComplexF64[-2 0;0 1] rtol=3e-15
            @test full.integrals.spectral_trace_cumulative[:,k]≈-cumsum(problem.grids.wᴱ)/(2π) rtol=3e-15
            @test full.integrals.negative_energy_integral_scaled[1,k]≈2energy_weight rtol=3e-15
            @test full.integrals.negative_energy_integral_scaled[5,k]≈(3sqrt(5.0)-5)*energy_weight rtol=4e-15
        end
        expected_sheet=2energy_weight*problem.physical.g_s*sum(problem.grids.wᵏ)/problem.scales.L₀_m^2
        @test full.integrals.negative_phase_space_weight_per_m2[1]≈expected_sheet rtol=3e-15
        full.integrals.spectral_integral[1,1,1]=123
        @test solution.scba.green.A[1,1,1,1]==-2
        # Keep every complex coherence, rather than replacing the integral by populations.
        for e in 1:ne,k in 1:nk
            solution.scba.green.A[e,k,:,:]=ComplexF64[2 im;-im 2]
        end
        coherent=positivity_diagnostics(solution;integrate=true)
        @test coherent.integrals.spectral_integral[1,:,:]≈energy_weight*ComplexF64[2 im;-im 2] rtol=3e-15
        @test coherent.integrals.spectral_integral_eigenvalues[1,:]≈energy_weight*[1.0,3.0] rtol=3e-15
        # A literal budget at one solver unit: negative is measured/countable even
        # below the declared backward-error allowance; the state is never clipped.
        solution.scba.green.A.=0
        solution.scba.green.A[:,:,1,1].=-1e-16
        tiny=positivity_diagnostics(solution;integrate=true)
        @test tiny.witnesses[1].backward_error==64*2*eps(Float64)
        @test tiny.witnesses[1].absolute_defect==1e-16 && tiny.witnesses[1].ratio==0
        @test tiny.strict_threshold_satisfied[1]
        @test tiny.integrals.negative_block_count[1]==ne*nk
        @test tiny.integrals.beyond_backward_error_block_count[1]==0
        @test solution.scba.green.A[1,1,1,1]==-1e-16
        for scale in (1e-6,1.0,1e6)
            solution.scba.green.A.=0
            solution.scba.green.A[:,:,1,1].=-2scale
            solution.scba.green.A[:,:,2,2].=scale
            scaled=positivity_diagnostics(solution)
            budget=64*2*eps(Float64)*max(2scale,1.0)
            @test scaled.witnesses[1].ratio≈(2scale-budget)/(2scale) rtol=3e-15
            @test scaled.strict_absolute_limits_scaled[1]≈budget+solution.options.tolerances.r_PSD*2scale rtol=3e-15
        end
        # Strict > tie selection has different historical scan orders; retain both.
        solution.scba.green.A.=0
        solution.scba.green.A[:,:,1,1].=1;solution.scba.green.A[:,:,2,2].=1
        solution.scba.green.A[2,1,1,1]=-2;solution.scba.green.A[1,2,1,1]=-2
        a=positivity_diagnostics(solution).witnesses[1]
        b=positivity_diagnostics(solution;integrate=true).witnesses[1]
        @test (a.energy_index,a.momentum_index)==(1,2)
        @test (b.energy_index,b.momentum_index)==(2,1)
        solution.scba.green.A[1,1,1,1]=NaN
        invalid=positivity_diagnostics(solution;integrate=true)
        @test invalid.witnesses[1].status===:error && isnan(invalid.witnesses[1].minimum_eigenvalue)
        @test !invalid.strict_threshold_satisfied[1]
        @test invalid.integrals.invalid_block_count[1]==1
        @test isnan(invalid.integrals.negative_energy_integral_scaled[1,1])
        @test !ismutabletype(typeof(invalid.witnesses[1]))
        @test invalid.witnesses[1].hermiticity_defect==0 # preserved legacy default, status=:error
        @test !positivity_diagnostics(solution).strict_threshold_satisfied[1]
    end
end
end
