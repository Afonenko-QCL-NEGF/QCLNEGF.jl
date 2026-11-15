module DirectStaticRetardedVariant
include("../support/common.jl")
const BN = QCLNEGF

function fixture()
    base = tutorial_numerics()
    changes = (N_z=24, N_b=2, P_basis=1, N_E=9, N_k=2, N_φ=4, N_qz=9)
    numerical = NumericalParameters(;
        (name => get(changes, name, getfield(base, name)) for
         name in fieldnames(NumericalParameters))...,
    )
    physical = reference_parameters(V_period=0u"mV", Tᴸ=70u"K",
        ΔV_alloy=0.60u"eV", Ω₀=4.5168e-29u"m^3")
    problem = build_problem(; physical, numerical,
        scattering=ScatteringOptions(alloy=true), localization=:none)
    ne, nk, nb = numerical.N_E, numerical.N_k, numerical.N_b
    GR = zeros(ComplexF64, ne, nk, nb, nb)
    # Independent rational inverse of H=[0 1+i; 1-i 2]. Its poles are 1±sqrt(3).
    for e=1:ne, m=1:nk
        z = problem.grids.ε[e] - .17m + (.2+.1m)*im
        GR[e,m,:,:] .= ComplexF64[z-2 1+im; 1-im z] / (z*z-2z-2)
    end
    A = spectral_function(GR)
    GL = 0.37im .* A
    GG = -0.63im .* A
    green = GreenState(GR, GL, GG, A, ones(ne,nk), ones(ne,nk))
    return problem, green
end

function replace_kernels(problem, K, names; qscale=1.0)
    tensors = Dict(name=>K./qscale for name in names)
    scales = Dict(name=>Float64(qscale) for name in names)
    kernels = KernelSet(tensors, scales, copy(problem.kernels.Fᴸᴼ), collect(names))
    return NEGFProblem(problem.physical, problem.numerical, problem.scattering,
        problem.scales, problem.grids, problem.profiles, problem.basis, kernels,
        problem.W₊ᴱᵖ, problem.W₋ᴱᵖ, problem.W₊ᴸᴼ, problem.W₋ᴸᴼ,
        problem.energy_shift_discretization, problem.models)
end

function production_options(mode; drop=false)
    return ProductionOptions(worker_count=1,
        algorithms=AlgorithmOptions(elastic_retarded=mode, hilbert=:direct,
            kernel_build=:direct, retarded_real_part=drop ? :drop : :kramers_kronig))
end

# No Core contraction helper or flattened tensor is used by this oracle.
function covariance_oracle(X, vertices, weights)
    ne,nk,nb,_ = size(X)
    result = zeros(ComplexF64, size(X))
    for e=1:ne, m=1:nk, a=1:nb, b=1:nb
        for mp=1:nk, V in vertices, c=1:nb, d=1:nb
            coefficient = (m+mp)/(nk+1)
            result[e,m,a,b] += coefficient*weights[mp]*V[a,c]*X[e,mp,c,d]*conj(V[b,d])
        end
    end
    return result
end

# Literal six-index oracle for actual assembled kernels, independent from BLAS.
function actual_kernel_oracle(problem, name, X)
    ne,nk,nb,_ = size(X)
    K = problem.kernels.K[name]
    out = zeros(ComplexF64, size(X))
    for e=1:ne, m=1:nk, a=1:nb, b=1:nb
        for mp=1:nk, c=1:nb, d=1:nb
            out[e,m,a,b] += problem.kernels.qᴷ[name]*problem.grids.wᵏ[mp]*
                           K[m,mp,a,c,d,b]*X[e,mp,c,d]
        end
    end
    return out
end

function broadening(R)
    out = similar(R)
    for e in axes(R,1), m in axes(R,2)
        block = Matrix(view(R,e,m,:,:))
        out[e,m,:,:] .= im .* (block-block')
    end
    return out
end

@testset "Explicit direct static SCBA closure" begin
    original, green = fixture()
    vertices = [ComplexF64[1 .2+.3im; .2-.3im .7],
                ComplexF64[.5 .4im; -.4im 1.2]]
    @test norm(vertices[1]*vertices[2]-vertices[2]*vertices[1]) > .1
    nk, nb = original.numerical.N_k, original.numerical.N_b
    K = zeros(ComplexF64,nk,nk,nb,nb,nb,nb)
    for m=1:nk, mp=1:nk, a=1:nb, c=1:nb, d=1:nb, b=1:nb
        K[m,mp,a,c,d,b] = (m+mp)/(nk+1)*
            sum(V[a,c]*conj(V[b,d]) for V in vertices)
    end
    names = (:acoustic,:impurity,:IFR,:alloy)
    problem = replace_kernels(original,K,names;qscale=7.25)
    po = production_options(:direct_green)
    cache = build_production_cache(problem; options=po)
    ref, rref = BN._scattering_candidate(problem,green;
        elastic_retarded=:direct_green,return_roundoff=true)
    prod, rprod = BN._scattering_candidate_production(problem,green,cache,po;
        return_roundoff=true)
    truthR = covariance_oracle(green.Gᴿ,vertices,problem.grids.wᵏ)
    truthL = covariance_oracle(green.Gˡ,vertices,problem.grids.wᵏ)
    truthG = covariance_oracle(green.Gᵍ,vertices,problem.grids.wᵏ)
    truthA = covariance_oracle(green.A,vertices,problem.grids.wᵏ)
    for name in names, candidate in (ref,prod)
        @test candidate[name].Σᴿ ≈ truthR rtol=5e-13 atol=2e-12
        @test candidate[name].Σˡ ≈ truthL rtol=5e-13 atol=2e-12
        @test candidate[name].Σᵍ ≈ truthG rtol=5e-13 atol=2e-12
        @test broadening(candidate[name].Σᴿ) ≈ truthA rtol=5e-13 atol=2e-12
        @test im.*(candidate[name].Σᵍ-candidate[name].Σˡ) ≈ truthA rtol=5e-13 atol=2e-12
        for e in axes(truthA,1), m in axes(truthA,2)
            @test eigmin(Hermitian(Matrix(view(truthA,e,m,:,:)))) >= -2e-12
        end
    end
    # No Hilbert reconstruction was applied in this all-static fixture.
    @test rref == rprod == 0.0
    @test BN._static_contraction(K,im.*green.Gᴿ,problem.grids.wᵏ,1.0) ≈ im.*truthR

    @testset "Reference final audit retains the selected map" begin
        algorithms = AlgorithmOptions(solver_backend=:educational,energy_shift=:dense,
            hilbert=:direct,contraction=:literal,kernel_build=:direct,
            elastic_retarded=:direct_green)
        options = SolverOptions(max_scba=1,max_poisson=1)
        embedding,plus,minus = embedding_self_energy(problem,green)
        state = SCBAResult(green,ref,embedding,plus,minus,SCBAIteration[],
            false,:fixture,:unresolved,BN._solver_restart_contract(options,algorithms))
        solution = NEGFSolution(problem,options,zeros(problem.numerical.N_z),
            BN._electron_density_bar(problem,green.Gˡ),state,OuterIteration[],
            Dict{Symbol,Any}(:elastic_retarded=>"direct_green"),
            ConvergenceReport(false,Dict{Symbol,Float64}(),String[]),false,:fixture)
        @test BN._algorithm_options_from_solution(solution).elastic_retarded === :direct_green
        metadata = BN.QCLNumerics._acceptance_metric_metadata(:r_roundoff,solution)
        @test metadata["applicability"]["applies"] === false
        audit = BN._audit_candidate_for_solution(solution)
        for name in names
            @test audit.scattering[name].Σᴿ == ref[name].Σᴿ
        end
        @test_throws ArgumentError solve_scba(problem,zeros(problem.numerical.N_z);
            options,initial=state)
        @test_throws ArgumentError solve_scba(problem,zeros(problem.numerical.N_z);
            options,elastic_retarded=:invalid)
    end

    @testset "Actual four elastic tensors preserve the retarded identity" begin
        actual = BN._scattering_candidate(original,green; elastic_retarded=:direct_green)
        op = production_options(:direct_green)
        actualcache = build_production_cache(original; options=op)
        fast = BN._scattering_candidate_production(original,green,actualcache,op)
        for name in names
            truth = actual_kernel_oracle(original,name,green.A)
            for candidate in (actual,fast)
                @test broadening(candidate[name].Σᴿ) ≈ truth rtol=2e-12 atol=2e-12
                @test im.*(candidate[name].Σᵍ-candidate[name].Σˡ) ≈ truth rtol=2e-12 atol=2e-12
            end
        end
    end

    @testset "Independent scalar semicircle through the direct route" begin
        s = .3
        Ks = zeros(ComplexF64,size(K))
        for m=1:nk,mp=1:nk,a=1:nb,b=1:nb
            Ks[m,mp,a,a,b,b] = s^2/sum(problem.grids.wᵏ)
        end
        scalar = replace_kernels(original,Ks,(:acoustic,))
        GR = zeros(ComplexF64,size(green.Gᴿ))
        roots = ComplexF64[]
        for (e,x) in enumerate(scalar.grids.ε)
            z = x + .15im
            g = 2/(z+sqrt(z-2s)*sqrt(z+2s))
            push!(roots,g)
            @test imag(g) < 0
            @test (z-s^2*g)*g ≈ 1 rtol=2e-13
            for m=1:nk,a=1:nb
                GR[e,m,a,a] = g
            end
        end
        A = spectral_function(GR)
        carrier = GreenState(GR,.37im.*A,-.63im.*A,A,
            ones(size(GR,1),nk),ones(size(GR,1),nk))
        candidate = BN._scattering_candidate(scalar,carrier; elastic_retarded=:direct_green)
        cp = production_options(:direct_green)
        cc = build_production_cache(scalar;options=cp)
        fast = BN._scattering_candidate_production(scalar,carrier,cc,cp)
        for family in (candidate[:acoustic],fast[:acoustic])
            for e in eachindex(roots),m=1:nk,a=1:nb
                @test family.Σᴿ[e,m,a,a] ≈ s^2*roots[e] rtol=2e-13
                @test -2imag(family.Σᴿ[e,m,a,a]) ≈ -2s^2*imag(roots[e]) rtol=2e-13
            end
        end
    end

    @testset "Default and nonstatic channel remain unchanged" begin
        implicit = BN._scattering_candidate(original,green)
        explicit = BN._scattering_candidate(original,green;elastic_retarded=:hilbert)
        direct = BN._scattering_candidate(original,green;elastic_retarded=:direct_green)
        ph = production_options(:hilbert)
        ch = build_production_cache(original;options=ph)
        productionH = BN._scattering_candidate_production(original,green,ch,ph)
        pd = production_options(:direct_green)
        cd = build_production_cache(original;options=pd)
        productionD = BN._scattering_candidate_production(original,green,cd,pd)
        for name in original.kernels.enabled, field in (:Σᴿ,:Σˡ,:Σᵍ)
            @test getfield(implicit[name],field) == getfield(explicit[name],field)
        end
        for field in (:Σᴿ,:Σˡ,:Σᵍ)
            @test getfield(direct[:LO],field) == getfield(implicit[:LO],field)
            @test getfield(productionD[:LO],field) == getfield(productionH[:LO],field)
        end
        dropH,dropD = production_options(:hilbert;drop=true),production_options(:direct_green;drop=true)
        cacheH,cacheD = build_production_cache(original;options=dropH),build_production_cache(original;options=dropD)
        familiesH = BN._scattering_candidate_production(original,green,cacheH,dropH)
        familiesD = BN._scattering_candidate_production(original,green,cacheD,dropD)
        for name in original.kernels.enabled,field in (:Σᴿ,:Σˡ,:Σᵍ)
            @test getfield(familiesH[name],field) == getfield(familiesD[name],field)
        end
        # Finite-window Hilbert and direct static closure are intentionally distinct.
        finite_window = BN._scattering_candidate(problem,green)
        @test norm(ref[:acoustic].Σᴿ-finite_window[:acoustic].Σᴿ) > 1e-6
    end
end

@testset "Variant identity, validation and legacy Hilbert restart" begin
    options = SolverOptions()
    defaults = AlgorithmOptions()
    direct = AlgorithmOptions(elastic_retarded=:direct_green,kernel_build=:direct)
    @test defaults.elastic_retarded === :hilbert
    @test algorithm_manifest(direct)["elastic_retarded"] == "direct_green"
    @test algorithm_impact(direct) === :controlled_numerical
    @test_throws ArgumentError BN._check_algorithm_options(AlgorithmOptions(elastic_retarded=:unknown))
    @test_throws ArgumentError BN._check_algorithm_options(AlgorithmOptions(
        elastic_retarded=:direct_green,contraction=:low_rank,low_rank_relative_tolerance=1e-4))
    @test_throws ArgumentError BN._check_algorithm_options(AlgorithmOptions(
        elastic_retarded=:direct_green,contraction=:low_rank,low_rank_maximum_rank=2))
    old = BN._solver_restart_contract(options,defaults)
    delete!(old["algorithms"],"elastic_retarded")
    original = deepcopy(old)
    @test BN._check_exact_restart_contract(old,options,defaults) === nothing
    @test old == original
    @test_throws ArgumentError BN._check_exact_restart_contract(old,options,direct)
    @test BN._check_exact_restart_contract(BN._solver_restart_contract(options,direct),options,direct) === nothing
end

end
