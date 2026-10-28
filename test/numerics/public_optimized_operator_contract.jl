module PublicOptimizedOperatorContractTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")

@testset "Public mix admits only pairwise independent writable target storage" begin
    fields=(:Σᴿ,:Σˡ,:Σᵍ)
    family(a,b,c)=SelfEnergyFamily((fill(ComplexF64(x),2,1,1,1) for x in (a,b,c))...)
    snapshot(f)=Tuple(copy(getfield(f,name)) for name in fields)
    unchanged(f,before)=all(isequal(getfield(f,name),before[i]) for (i,name) in enumerate(fields))
    options=ProductionOptions(worker_count=1)
    for overlap in ((1,2),(1,3),(2,3),(1,2,3))
        arrays=[fill(ComplexF64(x),2,1,1,1) for x in (2,4,8)]
        for index in overlap[2:end]
            arrays[index]=arrays[first(overlap)]
        end
        target=SelfEnergyFamily(arrays...);candidate=family(10,20,40)
        before=snapshot(target);candidate_before=snapshot(candidate)
        @test_throws ArgumentError production_mix_selfenergy!(target,candidate,0.5;options)
        @test unchanged(target,before)
        @test unchanged(candidate,candidate_before)
    end
    # Candidate components may share storage because the operator only reads them.
    target=family(2,4,8);shared=fill(ComplexF64(10),2,1,1,1)
    candidate=SelfEnergyFamily(shared,shared,shared);candidate_before=snapshot(candidate)
    @test production_mix_selfenergy!(target,candidate,0.5;options)===target
    @test all(==(6.0+0im),target.Σᴿ) && all(==(7.0+0im),target.Σˡ) && all(==(9.0+0im),target.Σᵍ)
    @test unchanged(candidate,candidate_before)
    healthy=family(2,4,8);candidate=family(10,20,40)
    @test production_mix_selfenergy!(healthy,candidate,0.5;options)===healthy
    @test (first(healthy.Σᴿ),first(healthy.Σˡ),first(healthy.Σᵍ))==(6+0im,12+0im,24+0im)
    # Adjacent admission: every cross-family storage pair and every component
    # shape is checked before any of the three target arrays is written.
    for left in fields,right in fields
        target=family(2,4,8);candidate=family(10,20,40)
        setfield!(candidate,right,getfield(target,left))
        before=snapshot(target);candidate_before=snapshot(candidate)
        @test_throws ArgumentError production_mix_selfenergy!(target,candidate,0.5;options)
        @test unchanged(target,before) && unchanged(candidate,candidate_before)
    end
    for side in (:target,:candidate),field in fields
        target=family(2,4,8);candidate=family(10,20,40)
        setfield!(side===:target ? target : candidate,field,zeros(ComplexF64,3,1,1,1))
        before=snapshot(target);candidate_before=snapshot(candidate)
        @test_throws DimensionMismatch production_mix_selfenergy!(target,candidate,0.5;options)
        @test unchanged(target,before) && unchanged(candidate,candidate_before)
    end
    for alpha in (0.0,-0.5,1.01,NaN,Inf)
        target=family(2,4,8);candidate=family(10,20,40)
        before=snapshot(target);candidate_before=snapshot(candidate)
        @test_throws ArgumentError production_mix_selfenergy!(target,candidate,alpha;options)
        @test unchanged(target,before) && unchanged(candidate,candidate_before)
    end
    target=family(2,4,8);candidate=family(10,20,40)
    @test production_mix_selfenergy!(target,candidate,1.0;options)===target
    @test snapshot(target)==snapshot(candidate)
end

@testset "Public LO collision keeps every kernel index, shift sign and Bose factor" begin
    available=isdefined(QCLNEGF,:production_lo_selfenergies)
    @test available
    if available
        base=public_state_problem()
        n=NumericalParameters(; (name=>(name===:N_b ? 2 : getfield(base.numerical,name))
            for name in fieldnames(NumericalParameters))...)
        scattering=ScatteringOptions(LO=true,acoustic=false,impurity=false,IFR=false,alloy=false)
        problem=build_problem(physical=base.physical,numerical=n,scattering=scattering,
            scales=base.scales,validate_static=false)
        shape=(n.N_E,n.N_k,n.N_b,n.N_b)
        GL=reshape(ComplexF64[0.1j+0.2im*j for j in 1:prod(shape)],shape)
        GG=reshape(ComplexF64[-0.15j-0.3im*j for j in 1:prod(shape)],shape)
        green=GreenState(zero(GL),GL,GG,zero(GL),ones(n.N_E,n.N_k),ones(n.N_E,n.N_k))
        saved=deepcopy(green)
        occupation=1/expm1(ustrip(u"eV",problem.physical.ħωᴸᴼ)/
            ustrip(u"eV",CODATA.kᴮₑᵥ*problem.physical.Tᴸᴼ))
        # Literal formula in (e,m,a,b; ep,mp,c,d) order. It deliberately calls no
        # shift/contraction/Bose helper from either reference or optimized Core.
        function literal_lo(G,plus,minus,population=occupation)
            out=zeros(ComplexF64,shape)
            K=problem.kernels.K[:LO]; q=problem.kernels.qᴷ[:LO]
            for e in 1:n.N_E,m in 1:n.N_k,a in 1:n.N_b,b in 1:n.N_b
                for ep in 1:n.N_E,mp in 1:n.N_k,c in 1:n.N_b,d in 1:n.N_b
                    out[e,m,a,b]+=q*problem.grids.wᵏ[mp]*K[m,mp,a,c,d,b]*
                        ((population+1)*plus[e,ep]+population*minus[e,ep])*G[ep,mp,c,d]
                end
            end
            out
        end
        expected_less=literal_lo(GL,problem.W₊ᴸᴼ,problem.W₋ᴸᴼ)
        expected_greater=literal_lo(GG,problem.W₋ᴸᴼ,problem.W₊ᴸᴼ)
        @test norm(expected_less)>0 && norm(expected_greater)>0
        for contraction in (:dense_blas,:literal)
            options=ProductionOptions(worker_count=1,algorithms=AlgorithmOptions(contraction=contraction))
            cache=build_production_cache(problem;options)
            result=production_lo_selfenergies(problem,green,cache;options)
            @test result.lesser≈expected_less rtol=3e-14
            @test result.greater≈expected_greater rtol=3e-14
            @test isequal(green.Gˡ,saved.Gˡ) && isequal(green.Gᵍ,saved.Gᵍ)
            result.lesser.=123
            @test isequal(green.Gˡ,saved.Gˡ)
            stale=public_state_problem()
            @test_throws ArgumentError production_lo_selfenergies(stale,green,cache;options)
            fixed_models=PhysicalModelOptions(lo_population=:fixed_occupation,lo_fixed_occupation=0.375)
            fixed=typeof(problem)((name===:models ? fixed_models : getfield(problem,name)
                for name in fieldnames(typeof(problem)))...)
            nonthermal=production_lo_selfenergies(fixed,green,cache;options)
            @test nonthermal.lesser≈literal_lo(GL,problem.W₊ᴸᴼ,problem.W₋ᴸᴼ,0.375) rtol=3e-14
            @test nonthermal.greater≈literal_lo(GG,problem.W₋ᴸᴼ,problem.W₊ᴸᴼ,0.375) rtol=3e-14
        end
    end
end

@testset "Public optimized resolvent and family mix match literal scalar and matrix truths" begin
    available=all(name->isdefined(QCLNEGF,name),(:production_retarded_green,
        :production_mix_selfenergy!,:fermi_occupation,:sampled_linewidth_diagnostics,
        :production_lo_selfenergies))
    @test available
    if available
        options=ProductionOptions(worker_count=1)
        energies=[-0.5,0.0,0.5]
        h=reshape(ComplexF64[0.125],1,1,1)
        sigma=fill(-0.25im,3,1,1,1)
        out=production_retarded_green(energies,h,sigma;options)
        @test out.G_R[:,1,1,1]≈[inv(e-0.125+0.25im) for e in energies] rtol=2e-15
        @test out.condition_number==ones(3,1)
        @test out.dyson_scale[:,1]≈abs.(energies .- 0.125 .+ 0.25im) rtol=2e-15
        @test_throws DimensionMismatch production_retarded_green(energies,h,sigma[1:2,:,:,:];options)
        # A noncommuting two-level literal adjugate is independent of the production inverse.
        hm=reshape(ComplexF64[0.1 0.02im;-0.02im -0.1],1,2,2)
        sm=reshape(ComplexF64[-0.03im 0;0 -0.07im],1,1,2,2)
        actual=production_retarded_green([0.05],hm,sm;options)
        a,d=-0.05+0.03im,0.15+0.07im
        b,c=-0.02im,0.02im
        expected=ComplexF64[d -b;-c a]/(a*d-b*c)
        @test actual.G_R[1,1,:,:]≈expected rtol=5e-15
        # The analytic causal root is supplied once, not iterated by a solver.
        z=0.2+0.3im; coupling=0.04
        root=(z-sqrt(z*z-4coupling))/(2coupling)
        imag(root)<0 || (root=(z+sqrt(z*z-4coupling))/(2coupling))
        fixed=production_retarded_green([real(z)],zeros(ComplexF64,1,1,1),
            fill(coupling*root,1,1,1,1);η=imag(z),options)
        @test only(fixed.G_R)≈root rtol=5e-15
        family(a,b,c)=SelfEnergyFamily(fill(ComplexF64(a),1,1,1,1),
            fill(ComplexF64(b),1,1,1,1),fill(ComplexF64(c),1,1,1,1))
        old=family(1,2im,-3im); candidate=family(5,6im,-7im)
        saved=deepcopy(candidate)
        @test production_mix_selfenergy!(old,candidate,0.25;options)===old
        @test (only(old.Σᴿ),only(old.Σˡ),only(old.Σᵍ))==(2.0+0im,3.0im,-4.0im)
        @test isequal(candidate.Σᴿ,saved.Σᴿ) && isequal(candidate.Σˡ,saved.Σˡ)
        @test production_selfenergy_residual(old,candidate;worker_count=1)≈0.6 rtol=2e-14
        @test_throws ArgumentError production_mix_selfenergy!(old,candidate,0.0;options)
        malformed=SelfEnergyFamily(zeros(ComplexF64,1,1,1,1),zeros(ComplexF64,2,1,1,1),
            zeros(ComplexF64,1,1,1,1))
        @test_throws DimensionMismatch production_mix_selfenergy!(old,malformed,0.5;options)
        @test_throws ArgumentError production_mix_selfenergy!(old,old,0.5;options)
    end
end

@testset "Positive thermal Fermi wrapper preserves its saturation and private zero-temperature semantics" begin
    available=isdefined(QCLNEGF,:fermi_occupation)
    @test available
    if available
        for (energy,want) in ((0.0,0.5),(log(3.0),0.25),(-log(3.0),0.75),(41.0,0.0),(-41.0,1.0))
            @test fermi_occupation(energy,0.0,1.0)≈want rtol=2e-15
        end
        @test fermi_occupation(40.0,0.0,1.0)>0
        for temperature in (0.0,-1.0,NaN,Inf)
            @test_throws ArgumentError fermi_occupation(0.0,0.0,temperature)
        end
        @test isnan(QCLNEGF.QCLNumerics._fermi(0.0,0.0,0.0))
        @test QCLNEGF.QCLNumerics._fermi(-1.0,0.0,0.0)==1.0
        @test QCLNEGF.QCLNumerics._fermi(1.0,0.0,0.0)==0.0
        bath=PhysicalModelExtensions.HotLOPhononBath(0.036u"eV",0.0u"K",1.0u"ps")
        @test PhysicalModelExtensions.equilibrium_lo_population(bath).occupation==0.0
    end
end

@testset "Sampled modal linewidth keeps Rayleigh bounds and unavailable status reasons" begin
    available=isdefined(QCLNEGF,:sampled_linewidth_diagnostics)
    @test available
    if available
        A=zeros(ComplexF64,3,1,2,2)
        A[:,:,1,1].=1; A[:,:,2,2].=2
        lesser=zeros(ComplexF64,size(A)); greater=copy(lesser)
        greater[:,:,1,1].=-2im; greater[:,:,2,2].=-4im
        total=SelfEnergyFamily(copy(lesser),lesser,greater)
        policy=SCBAPhysicsMarkerPolicy(max_spectral_blocks=3)
        result=sampled_linewidth_diagnostics(A,total,[-1.0,0.0,1.0],[0.5,1.0,0.5],[1.0];policy)
        @test result.status===:available
        @test (result.q10,result.q50,result.q90)==(2.0,4.0,4.0)
        @test result.underresolved==0 && result.blocks==3
        @test result.sampling_method===:spectral_midpoint_quantiles_A_eigenmodes_v1
        # Rayleigh bounds of diagonal gamma, independent of the marker implementation.
        @test 2.0<=result.q10<=result.q50<=result.q90<=4.0
        for (input,energy,reason) in ((zero(A),[-1.0,0.0,1.0],:empty_spectrum),
            (fill(ComplexF64(NaN),size(A)),[-1.0,0.0,1.0],:nonfinite_spectrum),
            (A,[0.0,0.0,0.0],:invalid_energy_grid),(-A,[-1.0,0.0,1.0],:negative_spectral_weight))
            unavailable=sampled_linewidth_diagnostics(input,total,energy,[0.5,1.0,0.5],[1.0];policy)
            @test unavailable.status===reason && unavailable.blocks==0
            @test all(isnan,(unavailable.q10,unavailable.q50,unavailable.q90,unavailable.underresolved))
        end
        indefinite=copy(A);indefinite[:,:,1,1].=-1
        @test sampled_linewidth_diagnostics(indefinite,total,[-1.,0.,1.],[0.5,1.,0.5],[1.];policy).status===:indefinite_spectrum
        noncausal=SelfEnergyFamily(copy(lesser),copy(lesser),-greater)
        @test sampled_linewidth_diagnostics(A,noncausal,[-1.,0.,1.],[0.5,1.,0.5],[1.];policy).status===:noncausal_broadening
        bad=deepcopy(total);bad.Σᵍ[1,1,1,1]=NaN
        failed=sampled_linewidth_diagnostics(A,bad,[-1.,0.,1.],[0.5,1.,0.5],[1.];policy)
        @test failed.status===:nonfinite_broadening && isnan(failed.q50)
        @test_throws DimensionMismatch sampled_linewidth_diagnostics(A,total,[-1.,0.,1.],[1.],[1.];policy)
    end
end
end
