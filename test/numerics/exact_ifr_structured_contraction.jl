module ExactIFRStructuredContraction
include("../support/common.jl")
using Random
const N = QCLNEGF.QCLNumerics

# Independent literal Kraus tensor, with complex phases and signed vertices.
function manufactured()
    Nk,Nb,J = 3,3,4
    u = ComplexF64[1 im 0.3+0.7im -0.2im; 0.2im 1 0.7-0.1im 1; 0.4 0.3im 1 -0.5+0.4im]
    d = [-0.7,1.1,-0.2,0.5]
    S = [1.0 0.2 0.3; 0.2 0.8 0.1; 0.3 0.1 1.2]
    D = [d[j]*u[:,j]*u[:,j]' for j in 1:J]
    K = zeros(ComplexF64,Nk,Nk,Nb,Nb,Nb,Nb)
    for m in 1:Nk, mp in 1:Nk, j in 1:J, a in 1:Nb, c in 1:Nb, dd in 1:Nb, b in 1:Nb
        K[m,mp,a,c,dd,b] += S[m,mp]*D[j][a,c]*conj(D[j][b,dd])
    end
    q = maximum(abs,K)
    L,R = zeros(ComplexF64,Nb^2,J),zeros(ComplexF64,J,Nb^2)
    for j in 1:J, b in 1:Nb, a in 1:Nb
        pair = a+Nb*(b-1)
        L[pair,j] = d[j]^2/q*u[a,j]*conj(u[b,j])
        R[j,pair] = conj(u[a,j])*u[b,j]
    end
    op = N._SeparableIFRProductionKernel(L,R,ComplexF64.(S),Nk,Nb)
    return op,K./q,q,D,S
end

function kraus_oracle(G,D,S,w)
    result = zeros(ComplexF64,size(G))
    for e in axes(G,1), m in axes(G,2), mp in axes(G,2), V in D
        result[e,m,:,:] .+= w[mp]*S[m,mp] .* (V*Matrix(view(G,e,mp,:,:))*V')
    end
    return result
end

@testset "Exact IFR complex matrix units, CP and batched apply" begin
    op,K,q,D,S = manufactured()
    w = [0.2,0.7,1.1]
    # Every input matrix unit, not only real/diagonal Hermitian cases.
    choi = [zeros(ComplexF64,9,9) for m in 1:3, mp in 1:3]
    for mp in 1:3, c in 1:3, d in 1:3
        G = zeros(ComplexF64,7,3,3,3)
        for e in 1:7
            G[e,mp,c,d] = sin(e/3)+im*cos(e/5)
        end
        expected = kraus_oracle(G,D,S,w)
        literal = QCLNEGF._static_contraction(K,G,w,q)
        @test literal ≈ expected rtol=2e-13 atol=2e-14
        for chunk in (1,3,16)
            actual = production_static_contraction(op,G,w,q;energy_chunk=chunk)
            @test actual ≈ expected rtol=2e-13 atol=2e-14
            if chunk==1
                for m in 1:3, b in 1:3, a in 1:3
                    choi[m,mp][(c-1)*3+a,(d-1)*3+b] = actual[1,m,a,b]/G[1,mp,c,d]
                end
            end
        end
    end
    for C in choi
        @test norm(C-C') <= 3e-13*max(norm(C),1)
        @test eigmin(Hermitian((C+C')/2)) >= -3e-13*max(norm(C),1)
    end
    Random.seed!(0x1331F)
    G = randn(ComplexF64,7,3,3,3)
    for e in 1:7,m in 1:3
        X = Matrix(view(G,e,m,:,:)); G[e,m,:,:] .= X*X'
    end
    expected = kraus_oracle(G,D,S,w)
    for backend in (:threads,:blas), workers in (1,min(2,Base.Threads.nthreads(:default)))
        actual = production_static_contraction(op,G,w,q;energy_chunk=3,
            parallel_backend=backend,worker_count=workers)
        @test actual ≈ expected rtol=2e-13 atol=2e-14
        for e in 1:7,m in 1:3
            B = Matrix(view(actual,e,m,:,:))
            @test norm(B-B') <= 3e-13*max(norm(B),1)
            @test eigmin(Hermitian((B+B')/2)) >= -3e-13*max(norm(B),1)
        end
    end
    # Arbitrary retarded real part, with independent positive broadening.
    H = randn(ComplexF64,size(G))
    for e in 1:7,m in 1:3
        P = Matrix(view(H,e,m,:,:)); H[e,m,:,:] .= (P+P')/2
    end
    GR = H-0.5im*G
    SR = production_static_contraction(op,GR,w,q;energy_chunk=3)
    broadening = im .* (SR .- permutedims(conj.(SR),(1,2,4,3)))
    @test broadening ≈ expected rtol=3e-13 atol=3e-14
    before = randn(ComplexF64,size(G)); output = copy(before)
    N._production_static_contraction!(output,op,GR,w,q,3;α=0.7-0.2im,β=0.3)
    @test output ≈ (0.7-0.2im)*kraus_oracle(GR,D,S,w)+0.3*before rtol=3e-13 atol=3e-14
    # Repeated use and a shorter final chunk must not retain old buffer data.
    second = production_static_contraction(op,2*GR,w,q;energy_chunk=2)
    @test second ≈ 2*SR rtol=3e-13 atol=3e-14
end

function canonical_problem(;Nb=3,physical=reference_parameters(),validate_static=true)
    base = tutorial_numerics()
    changes = (;N_E=7,N_k=3,N_b=Nb)
    n = NumericalParameters(;
        (name=>hasproperty(changes,name) ? getproperty(changes,name) : getfield(base,name)
         for name in fieldnames(NumericalParameters))...)
    options = ProductionOptions(algorithms=AlgorithmOptions(hilbert=:direct),worker_count=1)
    p = build_problem_production(;physical,numerical=n,production_options=options,validate_static,
        scattering=ScatteringOptions(LO=false,acoustic=false,IFR=true,impurity=false,alloy=false)).problem
    return p,options
end

@testset "Canonical IFR recipe, source preservation and conservative fallback" begin
    p,options = canonical_problem()
    raw = interface_roughness_kernel(p.physical,p.scales,p.grids,p.basis)
    source = p.kernels.qᴷ[:IFR] .* p.kernels.K[:IFR]
    # Independent reference finite-angle construction versus shared production recipe.
    @test source ≈ raw rtol=3e-13 atol=3e-14*maximum(abs,raw)
    original = copy(p.kernels.K[:IFR])
    cache = build_production_cache(p;options)
    op = cache.kernels[:IFR]
    generic = QCLNEGF._production_kernel(p.kernels.K[:IFR])
    @test N._kernel_operator_bytes(op) < N._kernel_operator_bytes(generic)
    @test p.kernels.K[:IFR] == original
    @test 2*length(p.physical.interfaces) <= p.numerical.N_b^2
    G = randn(ComplexF64,7,3,3,3)
    actual = production_static_contraction(op,G,p.grids.wᵏ,p.kernels.qᴷ[:IFR];energy_chunk=3)
    @test actual ≈ QCLNEGF._static_contraction(p.kernels.K[:IFR],G,p.grids.wᵏ,p.kernels.qᴷ[:IFR]) rtol=4e-13 atol=4e-14
    # Exercise the canonical constructor itself with a nontrivial complex
    # phase frame, rather than only manufacturing its private factor type.
    U = Diagonal(exp.(im .* [0.2,0.6,-0.7]))
    fields = Dict{Symbol,Any}(:Φ=>p.basis.Φ*U, :χ=>p.basis.χ*U)
    for name in (:H₀,:Z,:M⁻¹,:T₊,:T₋)
        fields[name] = U'*getfield(p.basis,name)*U
    end
    phase_basis = BasisData((get(fields,name,getfield(p.basis,name))
                            for name in fieldnames(BasisData))...)
    phase_raw = N._interface_roughness_kernel_separable(p.physical,p.scales,p.grids,phase_basis)
    independent = interface_roughness_kernel(p.physical,p.scales,p.grids,phase_basis)
    @test phase_raw ≈ independent rtol=3e-13 atol=3e-14*maximum(abs,independent)
    phase_kernels = N._normalize_production_kernels(Dict(:IFR=>phase_raw),p.kernels.Fᴸᴼ)
    phase_problem = NEGFProblem(p.physical,p.numerical,p.scattering,p.scales,p.grids,
        p.profiles,phase_basis,phase_kernels,p.W₊ᴱᵖ,p.W₋ᴱᵖ,p.W₊ᴸᴼ,p.W₋ᴸᴼ,
        p.energy_shift_discretization,p.models)
    phase_op = build_production_cache(phase_problem;options).kernels[:IFR]
    @test N._kernel_operator_bytes(phase_op) < N._kernel_operator_bytes(generic)
    rotated = similar(G); gauge_expected = similar(G); positive = similar(G)
    for e in 1:7,m in 1:3
        X = Matrix(view(G,e,m,:,:))
        rotated[e,m,:,:] .= U'*X*U
        gauge_expected[e,m,:,:] .= U'*Matrix(view(actual,e,m,:,:))*U
        positive[e,m,:,:] .= U'*(X*X')*U
    end
    phase_result = production_static_contraction(phase_op,rotated,p.grids.wᵏ,
        phase_kernels.qᴷ[:IFR];energy_chunk=3)
    phase_literal = QCLNEGF._static_contraction(phase_kernels.K[:IFR],rotated,
        p.grids.wᵏ,phase_kernels.qᴷ[:IFR])
    @test phase_result ≈ phase_literal rtol=4e-13 atol=4e-14
    @test phase_result ≈ gauge_expected rtol=5e-13 atol=5e-14
    phase_positive = production_static_contraction(phase_op,positive,p.grids.wᵏ,
        phase_kernels.qᴷ[:IFR];energy_chunk=3)
    for e in 1:7,m in 1:3
        B = Matrix(view(phase_positive,e,m,:,:))
        @test norm(B-B') <= 4e-13*max(norm(B),1)
        @test eigmin(Hermitian((B+B')/2)) >= -4e-13*max(norm(B),1)
    end
    # A caller-modified tensor is not silently projected onto a nearby recipe.
    p.kernels.K[:IFR][1,1,1,1,1,1] += 0.01
    modified = build_production_cache(p;options).kernels[:IFR]
    expected = QCLNEGF._static_contraction(p.kernels.K[:IFR],G,p.grids.wᵏ,p.kernels.qᴷ[:IFR])
    @test production_static_contraction(modified,G,p.grids.wᵏ,p.kernels.qᴷ[:IFR]) ≈ expected rtol=3e-13 atol=3e-14
    @test N._kernel_operator_bytes(modified) == N._kernel_operator_bytes(generic)
    # Even a single-ulp custom change must not be silently accepted by a
    # tolerance-based rank/separability test. Preserve the source entry exactly.
    p.kernels.K[:IFR] .= original
    changed = nextfloat(real(original[1,1,1,1,1,1]))
    p.kernels.K[:IFR][1,1,1,1,1,1] = changed
    one_ulp = build_production_cache(p;options).kernels[:IFR]
    @test p.kernels.K[:IFR][1,1,1,1,1,1] == changed
    @test N._kernel_operator_bytes(one_ulp) == N._kernel_operator_bytes(generic)
    @test production_static_contraction(one_ulp,G,p.grids.wᵏ,p.kernels.qᴷ[:IFR]) ≈
          QCLNEGF._static_contraction(p.kernels.K[:IFR],G,p.grids.wᵏ,p.kernels.qᴷ[:IFR]) rtol=3e-13 atol=3e-14
    # This tiny operator-container fixture does not claim basis adequacy.
    small,small_options = canonical_problem(Nb=1,validate_static=false)
    small_cache = build_production_cache(small;options=small_options)
    @test N._kernel_operator_bytes(small_cache.kernels[:IFR]) ==
          N._kernel_operator_bytes(QCLNEGF._production_kernel(small.kernels.K[:IFR]))
    # If the actual finite-angle factor rounds to momentum independence, keep
    # the already cheaper existing reduced map and its smaller memory bound.
    flat,flat_options = canonical_problem(physical=reference_parameters(Λᴵᶠᴿ=1e-8u"nm"))
    flat_cache = build_production_cache(flat;options=flat_options)
    @test N._kernel_operator_bytes(flat_cache.kernels[:IFR]) == flat.numerical.N_b^4*sizeof(ComplexF64)
end
end
