module NonneutralPeriodicPoissonAssessment

include("../support/common.jl")
# Accepted carrier correction: /tmp/qcl-oracle-carrier-design.md, SHA256
# 17c6dda56a2223bde9d20e98928e659d9e60fb7d6615b72ff94dcf190b829ea1.
# Manufactured functional inputs only; never the S1 all-five model.
const ATOL=2e-14
const RTOL=2e-13
const eSI=1.602176634e-19
const hbarSI=1.054571817e-34
const m0SI=9.1093837139e-31
const eps0SI=8.8541878128e-12
const E0eV=0.1
const E0J=E0eV*eSI
const L0m=1e-8
const mr=0.067
const J0ref=eSI*E0J/(hbarSI*L0m^2)
const lambda_ref=eSI^2/(eps0SI*E0J*L0m)
const hopping=hbarSI^2/(2m0SI*mr*E0J*L0m^2)

function manufactured_problem(Nz, Nb, localization)
    layers=Layer[Layer(d=(Nz*10.0)*u"nm", material=:manufactured_uniform,
        x_Al=0.0, Eᶜ=0.0u"eV", mᶻᵣ=mr, m_parallelᵣ=mr,
        εᵣ=1.0, doped=true)]
    physical=PhysicalParameters(layers, (Nz*1.0e16)*u"m^-2", 1.0,
        0.0u"V/m", 0.0u"m", 0.0u"eV", 70.0u"K", 70.0u"K",
        12.9, 10.89, 0.0367u"eV", 2.0e8u"m^-1", 2.0e8u"m^-1",
        1.0e-10u"m", 9.0e-9u"m", 7.0u"eV", 5317.0u"kg/m^3",
        5240.0u"m/s", nothing, nothing, 2, typeof(1.0u"m")[])
    numerical=NumericalParameters(N_z=Nz, N_b=Nb, P_basis=1,
        E_min=-0.1u"eV", E_max=0.1u"eV", N_E=3, M_E=0.05u"eV",
        k_max=1.0e8u"m^-1", N_k=2, N_φ=4, qz_max=1.0e8u"m^-1",
        N_qz=3, η_seed=0.001u"eV")
    scattering=ScatteringOptions(LO=false, acoustic=false, impurity=false,
        IFR=false, alloy=false)
    problem=build_problem(physical=physical, numerical=numerical,
        scattering=scattering, scales=ScaleSystem(E₀=0.1u"eV", L₀=10.0u"nm"),
        localization=localization, validate_static=true,
        physical_models=PhysicalModelOptions())
    @test validate_problem(problem).passed
    g,p,s=problem.grids,problem.profiles,problem.scales
    @test g.x ≈ collect(0.5:1.0:(Nz-0.5)) atol=ATOL rtol=RTOL
    @test g.wˣ ≈ ones(Nz) atol=ATOL rtol=RTOL
    @test g.ε == [-1.0,0.0,1.0]
    @test g.wᴱ == [0.5,1.0,0.5]
    @test g.trusted_energy == BitVector([false,true,false])
    @test g.κ ≈ [0.0,1.0] atol=ATOL rtol=RTOL
    @test g.wᵏ ≈ [1,3]./(16π) atol=ATOL rtol=RTOL
    @test g.qᶻ ≈ [-1.0,0.0,1.0] atol=ATOL rtol=RTOL
    @test g.wᑫᶻ ≈ [0.5,1.0,0.5] atol=ATOL rtol=RTOL
    @test g.φ ≈ [0,π/2,π,3π/2] atol=ATOL rtol=RTOL
    @test p.Eᶜ == zeros(Nz)
    @test p.mᶻᵣ == fill(mr,Nz) && p.m_parallelᵣ == fill(mr,Nz)
    @test p.εᵣ == ones(Nz) && p.x_Al == zeros(Nz)
    @test p.layer_index == ones(Int,Nz)
    @test p.Nᴰ ≈ ones(Nz) atol=ATOL rtol=RTOL
    @test isempty(problem.kernels.enabled) && isempty(problem.kernels.K) && isempty(problem.kernels.qᴷ)
    @test problem.kernels.Fᴸᴼ == zeros(ComplexF64,3,Nb,Nb)
    I3=Matrix{Float64}(I,3,3)
    @test problem.W₊ᴱᵖ == I3 && problem.W₋ᴱᵖ == I3
    delta=0.367
    @test problem.W₊ᴸᴼ ≈ [1-delta delta 0;0 1-delta delta;0 0 0] atol=ATOL rtol=RTOL
    @test problem.W₋ᴸᴼ ≈ [0 0 0;delta 1-delta 0;0 delta 1-delta] atol=ATOL rtol=RTOL
    @test problem.energy_shift_discretization == :nodal_linear
    defaults=PhysicalModelOptions()
    @test all(name->getfield(problem.models,name)==getfield(defaults,name),fieldnames(PhysicalModelOptions))
    @test s.E₀_eV/E0eV ≈ 1 atol=ATOL rtol=RTOL
    @test s.E₀_J/E0J ≈ 1 atol=ATOL rtol=RTOL
    @test s.L₀_m/L0m ≈ 1 atol=ATOL rtol=RTOL
    @test s.λ_P/lambda_ref ≈ 1 atol=ATOL rtol=RTOL
    @test ustrip(u"A/m^2",s.J₀)/J0ref ≈ 1 atol=ATOL rtol=RTOL
    return problem
end

function checked_plain_fields(arrays, condition, scale, shape)
    @assert all(X->X isa Array{ComplexF64,4} && size(X)==shape && all(isfinite,X), arrays)
    @assert size(condition)==size(scale)==shape[1:2]
    @assert all(x->isfinite(x)&&x>0,condition) && all(x->isfinite(x)&&x>0,scale)
    fields=(arrays...,condition,scale)
    @assert all(!Base.mightalias(fields[i],fields[j]) for i in eachindex(fields) for j in (i+1):length(fields))
end

function real_space_basis_inventory(problem)
    b=problem.basis
    I4=Matrix{ComplexF64}(I,4,4)
    H=hopping.*ComplexF64[2 -1 0 0;-1 2 -1 0;0 -1 2 -1;0 0 -1 2]
    seam=zeros(ComplexF64,4,4);seam[4,1]=-hopping
    @test b.localization==:real_space && b.Φ==I4
    @test b.χ ≈ I4 atol=ATOL rtol=RTOL
    @test b.centres ≈ [0.5,1.5,2.5,3.5] atol=ATOL rtol=RTOL
    @test b.spreads==zeros(4) && b.overlap_eigenvalues==ones(4) && b.κ_overlap==1
    @test b.H₀ ≈ H atol=ATOL rtol=RTOL
    @test b.Z ≈ Matrix(Diagonal([0.5,1.5,2.5,3.5])) atol=ATOL rtol=RTOL
    @test b.M⁻¹ ≈ I4./mr atol=ATOL rtol=RTOL
    @test b.T₊ ≈ seam atol=ATOL rtol=RTOL
    @test b.T₋ ≈ seam' atol=ATOL rtol=RTOL
    @test b.window_eigenvalues ≈ [2hopping*(1-cos(j*π/5)) for j in 1:4] atol=ATOL rtol=RTOL
    @test b.Egrid₀ ≈ 2hopping atol=ATOL rtol=RTOL
    @test b.r_orth==b.r_eigen==b.r_translation==b.r_translation_nearest==b.r_translation_two_pairs==0
    return H
end

function synthetic_routing_solution(problem,H)
    # Every metric sent below is a synthetic branch sentinel except the three
    # direct Poisson measurements. This envelope is never persisted/validated
    # as a scientific solution. Its spectral sum defect is explicitly nonzero.
    R=2π^2
    shape=(3,2,4,4)
    I4=Matrix{ComplexF64}(I,4,4)
    GR,Gl,Gg,A,SR,Sl,Sg=(zeros(ComplexF64,shape) for _ in 1:7)
    condition=ones(3,2);scale=fill(1/R,3,2)
    for e in 1:3,m in 1:2
        hm=H+hopping*problem.grids.κ[m]^2 .* I4
        GR[e,m,:,:].=-im*R.*I4
        Gl[e,m,:,:].=im*R.*I4
        Gg[e,m,:,:].=-im*R.*I4
        A[e,m,:,:].=2R.*I4
        SR[e,m,:,:].=problem.grids.ε[e].*I4-hm-im/R.*I4
        Sl[e,m,:,:].=im/R.*I4
        Sg[e,m,:,:].=-im/R.*I4
        D=problem.grids.ε[e].*I4-hm-SR[e,m,:,:]
        @test D ≈ im/R.*I4 atol=ATOL rtol=RTOL
        @test D*GR[e,m,:,:] ≈ I4 atol=ATOL rtol=RTOL
        @test GR[e,m,:,:]*Sl[e,m,:,:]*GR[e,m,:,:]' ≈ Gl[e,m,:,:] atol=ATOL rtol=RTOL
        @test GR[e,m,:,:]-GR[e,m,:,:]' ≈ Gg[e,m,:,:]-Gl[e,m,:,:] atol=ATOL rtol=RTOL
        @test A[e,m,:,:] ≈ im.*(GR[e,m,:,:]-GR[e,m,:,:]') atol=ATOL rtol=RTOL
        @test SR[e,m,:,:]-SR[e,m,:,:]' ≈ Sg[e,m,:,:]-Sl[e,m,:,:] atol=ATOL rtol=RTOL
        @test condition[e,m] ≈ cond(D) atol=ATOL rtol=RTOL
        @test scale[e,m] ≈ opnorm(D,Inf) atol=ATOL rtol=RTOL
        for X in (-im.*Gl[e,m,:,:],im.*Gg[e,m,:,:],A[e,m,:,:],
                  -im.*Sl[e,m,:,:],im.*Sg[e,m,:,:])
            @test X ≈ X' atol=ATOL rtol=RTOL
            @test eigmin(Hermitian(X)) >= -ATOL
        end
    end
    checked_plain_fields((GR,Gl,Gg,A,SR,Sl,Sg),condition,scale,shape)
    green=GreenState(GR,Gl,Gg,A,condition,scale)
    total=SelfEnergyFamily(SR,Sl,Sg)
    plus=SelfEnergyFamily(SR./2,Sl./2,Sg./2)
    minus=SelfEnergyFamily(SR./2,Sl./2,Sg./2)
    for name in (:Σᴿ,:Σˡ,:Σᵍ)
        @test getfield(plus,name)+getfield(minus,name)==getfield(total,name)
    end
    families=(total,plus,minus)
    @assert all(!Base.mightalias(getfield(families[i],name),getfield(families[j],other))
        for i in 1:3 for j in (i+1):3 for name in (:Σᴿ,:Σˡ,:Σᵍ) for other in (:Σᴿ,:Σˡ,:Σᵍ))
    # Independent weighted density=gs/(2pi)*2*(1/(4pi))*R=1.
    @test R/(2π^2)==1
    @test QCLNEGF._electron_density_bar(problem,Gl) ≈ ones(4) atol=ATOL rtol=RTOL
    @test QCLNEGF._spectral_sum_residual(problem,green) ≈ 4π-1 atol=ATOL rtol=RTOL
    scba=SCBAResult(green,Dict{Symbol,SelfEnergyFamily}(),total,plus,minus,
        SCBAIteration[],true,:converged,:strictly_converged)
    return NEGFSolution(problem,SolverOptions(),zeros(4),ones(4),scba,
        OuterIteration[],Dict{Symbol,Any}(:fixture_kind=>"synthetic_acceptance_routing"),
        ConvergenceReport(false,Dict{Symbol,Float64}(),String[]),false,:manufactured_routing)
end

# Frozen explicit base-gate inventory: never obtain sentinels from the target
# _solution_validation_limits, or a dropped production gate could hide itself.
const SYNTHETIC_GATE_KEYS=(
    :r_D,:r_A,:r_K,:r_Σ,:r_λ,:r_roundoff,:r_Jchange_scba,:r_population_scba,
    :r_PSD,:r_caus,:r_caus_sc,:r_caus_embedding,
    :r_selfenergy_identity,:r_selfenergy_identity_sc,:r_selfenergy_identity_embedding,
    :r_antihermitian_G,:r_antihermitian_Σ,:r_antihermitian_Σsc,:r_antihermitian_Σembedding,
    :r_sum,:r_J,:r_ImJ,:r_roundJ,:r_power,:r_ImP,:r_roundP,
    :r_P,:r_U,:r_n,:r_neutral,:r_ζ,:r_Jchange,:r_population,
    :r_tail_low,:r_tail_high,:r_tail_k,:r_edge_gamma,:r_edge_spectral)

@testset "Constant incompatible periodic source: direct diagnostics and synthetic routing" begin
    problem=manufactured_problem(4,4,:real_space)
    # Near-zero rP floor makes an approximate donor precondition inadequate.
    # Mandatory reviewer condition: literal ND==1 EXACT or stop/replan.
    exact_donors=problem.profiles.Nᴰ==ones(4)
    @test exact_donors
    exact_donors || error("DEV04 STOP/replan: donor profile is not exactly ones(4); retain literal RHS and tolerances")
    H=real_space_basis_inventory(problem)
    healthy=solve_periodic_poisson(problem,ones(4))
    incompatible=solve_periodic_poisson(problem,2 .* ones(4))
    U,zeta,rP,neutrality=healthy
    t=SolverTolerances()
    @test U ≈ zeros(4) atol=ATOL rtol=RTOL
    @test zeta==rP==neutrality==0
    @test abs(zeta)<=t.r_ζ && rP<=t.r_P && neutrality<=t.r_neutral
    badU,badzeta,badrP,badneutrality=incompatible
    @test dot(problem.grids.wˣ,badU) ≈ 0 atol=ATOL rtol=RTOL
    @test badU ≈ zeros(4) atol=ATOL rtol=RTOL
    @test abs(badzeta-lambda_ref/2)<=1e-10
    @test badrP ≈ 2/(2+1e-14) atol=ATOL rtol=RTOL
    @test badneutrality==1
    @test badrP>1/2 && badrP>t.r_P && badneutrality>t.r_neutral && abs(badzeta)>t.r_ζ
    # Forced mean-charge removal would feed the already measured healthy RHS.
    # Keep the incompatible expected diagnostics; they distinguish this variant.
    @test !isapprox(zeta,lambda_ref/2;atol=1e-10,rtol=0)
    @test !isapprox(rP,2/(2+1e-14);atol=ATOL,rtol=RTOL)
    @test neutrality!=1

    @testset "Synthetic routing inputs are not scientific acceptance" begin
        solution=synthetic_routing_solution(problem,H)
        sentinels=Dict{Symbol,Float64}(key=>0.0 for key in SYNTHETIC_GATE_KEYS)
        sentinels[:r_P]=rP;sentinels[:r_neutral]=neutrality;sentinels[:r_ζ]=abs(zeta)
        # r_U=0: the direct compatible U and published manufactured U are zero.
        @test solution.Uᴴ==U
        report=QCLNEGF._validate_solution_metrics(solution,sentinels)
        @test report.passed && isempty(report.messages)
        bad=copy(sentinels)
        measured=Dict(:r_P=>badrP,:r_neutral=>badneutrality,:r_ζ=>abs(badzeta))
        merge!(bad,measured)
        rejected=QCLNEGF._validate_solution_metrics(solution,bad)
        @test !rejected.passed
        failed=Set(first(split(message,'=')) for message in rejected.messages)
        @test failed==Set(["r_P","r_neutral","r_ζ"])
        @test length(rejected.messages)==3
        for key in (:r_P,:r_neutral,:r_ζ)
            isolated=copy(sentinels);isolated[key]=measured[key]
            report=QCLNEGF._validate_solution_metrics(solution,isolated)
            @test !report.passed
            @test length(report.messages)==1
            @test startswith(only(report.messages),string(key)*"=")
        end
        @test !solution.converged && !solution.report.passed
        @test solution.status==:manufactured_routing
    end
end
end
