module FinalStateAlignmentOracle

include("../support/common.jl")
# Complete manufactured FUNCTIONAL final-state oracle, not S1 or scientific
# acceptance. Core source: f5b73369a96d1225bd96fb34f5262c8dcf965f86.
# D1 failed exact-spacing premise is retained in /tmp/qcl-dev05-d1-full.jl.
# Approved D2 SHA256 eb9e1a8c39ad08bafa480970392d844c2c9fffdef02d847397ef87c52111bfc9.
# Independent D2 review SHA256 53c3e9b4061a6cdb42239d6d42f5b3f7a09b94eb36315b1ec22bfea103c4f55f.
# Breaks caught: stale stored metrics, predicted U substituted for published U,
# old G accepted at final U, or a required final gate dropped.
# Existing fresh_candidate_residual.jl remains the narrow-pole primitive;
# neither that primitive nor causal SCBA/maxiter tests are complete assembly.
const ATOL=2e-14 # approved, unchanged static-carrier guard only
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

# Literal inventory: do not derive from production limits (a removed gate
# must fail this test). There are 38 base and 3 supplied-candidate gates.
const GATES=(
    :r_D,:r_A,:r_K,:r_Σ,:r_λ,:r_roundoff,:r_Jchange_scba,:r_population_scba,
    :r_PSD,:r_caus,:r_caus_sc,:r_caus_embedding,:r_selfenergy_identity,
    :r_selfenergy_identity_sc,:r_selfenergy_identity_embedding,
    :r_antihermitian_G,:r_antihermitian_Σ,:r_antihermitian_Σsc,
    :r_antihermitian_Σembedding,:r_sum,:r_J,:r_ImJ,:r_roundJ,:r_power,
    :r_ImP,:r_roundP,:r_P,:r_U,:r_n,:r_neutral,:r_ζ,:r_Jchange,
    :r_population,:r_tail_low,:r_tail_high,:r_tail_k,:r_edge_gamma,
    :r_edge_spectral,:r_D_candidate,:r_A_candidate,:r_lambda_candidate)
const CANDIDATE_GATES=(:r_D_candidate,:r_A_candidate,:r_lambda_candidate)

# Operation ledger (six scalar blocks only): constants <=16 elementary real
# operations; source units use 5 Float64 operations (nm/ell/L0/Lp/h),
# stencil uses2 (h²,reciprocal), v64 uses2 (sqrt,reciprocal); all are guarded
# exact binary64 INPUT DATA. Outward ideal h factors add2 for tZ (square,
# divide),1 for U/u,1 for each n/ND,2 for chi (sqrt,divide),1 for L_h.
# Donor target/normalization uses<=8 real operations inside retained gamma16;
# density integrates six entries and SI conversion inside retained gamma32.
# active root <=16 plus cancellation factor (4c+z²)/|4c-z²|;
# inactive rational root <=16, with positive denominator guard. Every oracle
# operation below rounds OUTWARD at 384 bits. A Float64 literal is subsequently
# an exact input; its independent polynomial defect is retained, not assumed0.
# Float64 target comparison: gamma16 per complex triple product (<=12 rounded
# real operations), gamma32 per six-entry norm/weighted charge/projection,
# gamma64 per Dyson numerator+denominator+quotient. These are operation bounds,
# never physical tolerances or a bound for a general BLAS/FFT problem.
struct Ball
    lo::BigFloat
    hi::BigFloat
    function Ball(lo::BigFloat,hi::BigFloat)
        @assert lo<=hi
        new(lo,hi)
    end
end
rd(f)=setrounding(f,BigFloat,RoundDown)
ru(f)=setrounding(f,BigFloat,RoundUp)
B(x::Real)=Ball(BigFloat(x),BigFloat(x)) # exact for Float64/integers
B(x::AbstractString)=Ball(rd(()->parse(BigFloat,x)),ru(()->parse(BigFloat,x)))
Base.:+(a::Ball,b::Ball)=Ball(rd(()->a.lo+b.lo),ru(()->a.hi+b.hi))
Base.:-(a::Ball,b::Ball)=Ball(rd(()->a.lo-b.hi),ru(()->a.hi-b.lo))
Base.:-(a::Ball)=Ball(-a.hi,-a.lo)
function Base.:*(a::Ball,b::Ball)
    lows=[rd(()->x*y) for x in (a.lo,a.hi), y in (b.lo,b.hi)]
    highs=[ru(()->x*y) for x in (a.lo,a.hi), y in (b.lo,b.hi)]
    Ball(minimum(lows),maximum(highs))
end
function Base.:/(a::Ball,b::Ball)
    @assert b.lo>0 || b.hi<0
    a*Ball(rd(()->inv(b.hi)),ru(()->inv(b.lo)))
end
Base.sqrt(a::Ball)=(@assert a.lo>=0; Ball(rd(()->sqrt(a.lo)),ru(()->sqrt(a.hi))))
function Base.abs(a::Ball)
    a.lo>=0 && return a
    a.hi<=0 && return -a
    Ball(BigFloat(0),max(-a.lo,a.hi))
end
sq(a::Ball)=a*a
mid(a::Ball)=(a.lo+a.hi)/2
upper(a::Ball)=a.hi
lower(a::Ball)=a.lo
contains(a::Ball,x::Real)=a.lo<=BigFloat(x)<=a.hi
hull(a::Ball,b::Ball)=Ball(min(a.lo,b.lo),max(a.hi,b.hi))
expand(a::Ball,r::Ball)=Ball(max(BigFloat(0),rd(()->a.lo-r.hi)),ru(()->a.hi+r.hi))
gamma(n)=B(n)*B(eps(Float64)/2)/(B(1)-B(n)*B(eps(Float64)/2))
struct CBall
    re::Ball
    im::Ball
end
C(z::Number)=CBall(B(real(z)),B(imag(z)))
C(x::Ball)=CBall(x,B(0))
Base.:*(a::Ball,b::CBall)=C(a)*b
Base.:*(a::CBall,b::Ball)=a*C(b)
Base.:-(a::CBall)=CBall(-a.re,-a.im)
Base.:+(a::CBall,b::CBall)=CBall(a.re+b.re,a.im+b.im)
Base.:-(a::CBall,b::CBall)=CBall(a.re-b.re,a.im-b.im)
Base.:*(a::CBall,b::CBall)=CBall(a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re)
Base.conj(a::CBall)=CBall(a.re,-a.im)
Base.abs(a::CBall)=sqrt(sq(abs(a.re))+sq(abs(a.im)))
normbox(xs)=sqrt(foldl(+, (sq(abs(C(x))) for x in xs);init=B(0)))
function ratio(n::Ball,bn::Ball,d::Ball,bd::Ball)
    @assert lower(d)>upper(bd)>0
    Ball(max(BigFloat(0),rd(()->(n.lo-bn.hi)/(d.hi+bd.hi))),
         ru(()->(n.hi+bn.hi)/(d.lo-bd.hi)))
end
function check_enclosed(value, enclosure::Ball)
    @test isfinite(value)
    @test contains(enclosure,value)
end
function maxbox(xs)
    Ball(maximum(x.lo for x in xs),maximum(x.hi for x in xs))
end

function manufactured_problem(Nz, Nb, localization, a)
    layers=Layer[Layer(d=(Nz*10.0)*u"nm", material=:manufactured_uniform,
        x_Al=0.0, Eᶜ=0.0u"eV", mᶻᵣ=mr, m_parallelᵣ=mr,
        εᵣ=1.0, doped=true)]
    physical=PhysicalParameters(layers, (1.0e16/64)*u"m^-2", 1.0,
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
    @assert Nz==3 && Nb==1
    @assert unit(problem.physical.layers[1].d)==u"m" && unit(s.L₀)==u"m"
    @assert ustrip(u"m",problem.physical.layers[1].d)==a.units.ell
    @assert s.L₀_m==a.units.L0
    @assert length(problem.physical.layers)==1
    @assert ustrip(u"m",problem.physical.layers[1].d)/s.L₀_m==a.units.Lp
    @assert g.x==[(j-0.5)*a.units.h for j in 1:3]
    @assert g.wˣ==fill(a.units.h,3)
    @assert p.εᵣ==[1.0,1.0,1.0]
    # Input period Lp and rounded grid length3h are distinct; no equality
    # between them is asserted or used as a density/history normalizer.
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
    @test p.Nᴰ ≈ fill(Float64(mid(a.D/(B(3)*a.h))),Nz) atol=ATOL rtol=RTOL
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

function analytic_basis_phase(problem,a)
    Nb=problem.numerical.N_b
    full=ComplexF64[1/2 1/sqrt(2);1/sqrt(2) 0;1/2 -1/sqrt(2)]
    literal=full[:,1:Nb]
    overlaps=[dot(literal[:,a],problem.basis.Φ[:,a]) for a in 1:Nb]
    @assert all(z->isapprox(abs(z),1;atol=ATOL,rtol=RTOL),overlaps)
    phase=Diagonal(overlaps./abs.(overlaps))
    b=problem.basis
    @assert isapprox(b.Φ,literal*phase;atol=ATOL,rtol=RTOL)
    @test b.χ ≈ (literal./Float64(mid(sqrt(a.h))))*phase atol=ATOL rtol=RTOL
    @test b.localization == :none
    @test b.centres ≈ fill(Float64(mid(B(3)*a.h/B(2))),Nb) atol=ATOL rtol=RTOL
    @test b.spreads ≈ Float64(mid(a.h)).*[1/sqrt(2),1][1:Nb] atol=ATOL rtol=RTOL
    @test b.overlap_eigenvalues ≈ ones(Nb) atol=ATOL rtol=RTOL
    @test b.κ_overlap ≈ 1 atol=ATOL rtol=RTOL
    hoppingZ=Float64(mid(a.tZ))
    energies=[(2-sqrt(2))*hoppingZ,2hoppingZ][1:Nb]
    @test b.H₀ ≈ Matrix(Diagonal(ComplexF64.(energies))) atol=ATOL rtol=RTOL
    @test b.Z ≈ Float64(mid(a.h)).*(phase'*[1.5 -1/sqrt(2);-1/sqrt(2) 1.5][1:Nb,1:Nb]*phase) atol=ATOL rtol=RTOL
    @test b.M⁻¹ ≈ Matrix{ComplexF64}(I,Nb,Nb)./mr atol=ATOL rtol=RTOL
    seam=-hoppingZ.*[1/4 1/(2sqrt(2));-1/(2sqrt(2)) -1/2][1:Nb,1:Nb]
    @test b.T₊ ≈ phase'*seam*phase atol=ATOL rtol=RTOL
    @test b.T₋ ≈ b.T₊' atol=ATOL rtol=RTOL
    @test b.window_eigenvalues ≈ energies atol=ATOL rtol=RTOL
    @test b.Egrid₀ ≈ 2hoppingZ atol=ATOL rtol=RTOL
    @test isfinite(b.r_orth) && isfinite(b.r_eigen)
    @test b.r_translation == b.r_translation_nearest == b.r_translation_two_pairs == 0
    return phase
end

# All literals and interval comparisons are prepared before any fresh/final API.
function source_unit_literals()
    # RN source chain from original30nm/10nm arguments, never target grids.
    nm=1.0/1000000000.0
    ell=30.0*nm
    L0=10.0*nm
    Lp=ell/L0
    h=Lp/3.0
    h2=h*h
    qL=1.0/h2
    v=1.0/sqrt(3.0)
    @assert nm==0x1.12e0be826d695p-30 && ell==0x1.01b2b29a4692cp-25
    @assert L0==0x1.5798ee2308c3ap-27 && Lp==0x1.8000000000001p+1
    @assert h==0x1.0000000000001p+0
    @assert h2==0x1.0000000000002p+0 && qL==0x1.ffffffffffffcp-1
    @assert v==0x1.279a74590331dp-1
    @assert (1.0e16/64)*(L0*L0)==1/64
    @assert isfinite(h) && h>0 && isfinite(qL) && qL>0
    # At this particular source-derived qL the2/4 products and all row sums
    # are exact power-of-two scalings. Freeze before any Poisson call.
    @assert B(2qL).lo==(B(2)*B(qL)).lo
    @assert B(4qL).lo==(B(4)*B(qL)).lo
    @assert sum(abs,[2qL,-qL,-qL])==4qL
    (nm=nm,ell=ell,L0=L0,Lp=Lp,h=h,h2=h2,qL=qL,v=v)
end

function independent_literals()
    units=source_unit_literals()
    h=B(units.h) # exact real binary64 data, not rounded decimal1
    e,hbar,m0,eps0=B("1.602176634e-19"),B("1.054571817e-34"),B("9.1093837139e-31"),B("8.8541878128e-12")
    E0,L0,mass=B("0.1")*e,B("1e-8"),B("0.067")
    t0=sq(hbar)/(B(2)*m0*mass*E0*sq(L0))
    tZ=t0/sq(h)
    la=sq(e)/(eps0*E0*L0)
    tau=-tZ/B(4)
    c=sq(tZ)/B(8)
    piB=Ball(rd(()->BigFloat(π)),ru(()->BigFloat(π)))
    D=B(1)/B(64)
    H0=(B(2)-sqrt(B(2)))*tZ
    u=la*D*h/B(72)
    Uideal=[la*D*h/B(36)*B(j) for j in (-1,2,-1)]
    U=Float64.(mid.(Uideal))
    nideal=[D*B(j)/(B(4)*h) for j in (1,2,1)]
    n=Float64.(mid.(nideal))
    L=[2 -1 -1;-1 2 -1;-1 -1 2]
    @assert all(isfinite,U) && maximum(abs,U)>0
    @assert all(isfinite,n) && minimum(n)>0
    # Preserve gauge check with the independently bounded rounding of U_f;
    # exact ideal sum0 is not silently imposed on rounded published inputs.
    Uround=foldl(+, (abs(B(U[j])-Uideal[j]) for j in 1:3);init=B(0))
    Usum=foldl(+,B.(U);init=B(0))
    @assert upper(abs(Usum))<=upper(Uround)
    for j in 1:3
        idealRow=foldl(+, (B(L[j,k])*Uideal[k]/sq(h) for k in 1:3);init=B(0))
        idealRhs=la*(nideal[j]-D/(B(3)*h))
        @test contains(idealRow-idealRhs,0)
    end
    @test contains(foldl(+, (h*x for x in nideal);init=B(0))-D,0)
    @test length(GATES)==length(Set(GATES))==41
    # Rounded U/n are exact subsequent inputs; do not assume their ideal sums.
    (t0=t0,tZ=tZ,la=la,tau=tau,c=c,pi=piB,D=D,H0=H0,u=u,U=U,n=n,
     Uideal=Uideal,nideal=nideal,L=L,h=h,units=units,B=B(16)*sq(piB)*D)
end

function literal_state(problem,a,stale::Bool)
    GR,Gl,Gg,A=(Array{ComplexF64}(undef,3,2,1,1) for _ in 1:4)
    scales=Matrix{Float64}(undef,3,2)
    roots=Matrix{CBall}(undef,3,2)
    for e in 1:3,m in 1:2
        h=a.H0+(stale ? B(0) : a.u)+a.t0*B((m-1)^2)
        z=B(e-2)-h
        rad=B(4)*a.c-sq(z)
        active=(e==2 && m==1)
        # Fixed branch inventory, not chosen from measured target outputs.
        @assert active ? lower(rad)>0 : upper(rad)<0
        @assert lower(abs(z))>0
        cancellation=(B(4)*a.c+sq(z))/abs(rad)
        @assert isfinite(upper(cancellation))
        if active
            root=CBall(z/(B(2)*a.c),-sqrt(rad)/(B(2)*a.c))
        else
            denom=z+B(e==1 || e==2 ? -1 : 1)*sqrt(-rad)
            @assert lower(abs(denom))>0
            root=C(B(2)/denom)
        end
        roots[e,m]=root
        g=ComplexF64(Float64(mid(root.re)),Float64(mid(root.im)))
        @assert isfinite(g) && !iszero(g)
        # Correctly rounded literal enclosure; roots themselves never come
        # from target Dyson or SCBA routines.
        @test BigFloat(prevfloat(real(g)))<=root.re.lo<=root.re.hi<=BigFloat(nextfloat(real(g)))
        @test BigFloat(prevfloat(imag(g)))<=root.im.lo<=root.im.hi<=BigFloat(nextfloat(imag(g)))
        GR[e,m,1,1]=g
        av=-2imag(g)
        bv=active ? Float64(mid(a.B)) : 0.0 # spectral support proves these zeros
        @assert 0<=bv<=av
        A[e,m,1,1]=av
        Gl[e,m,1,1]=im*bv
        Gg[e,m,1,1]=-im*(av-bv)
        self=a.c*C(g)
        scales[e,m]=Float64(mid(abs(C(z)-self)))
        @assert scales[e,m]>0 && isfinite(scales[e,m])
        # Polynomial defect of actual Float64 literal, in outward >=256-bit
        # arithmetic. It is a construction error, never set to exact zero.
        defect=abs(C(a.c)*C(g)*C(g)-C(z)*C(g)+C(1))
        polybudget=gamma(16)*(abs(C(a.c)*C(g)*C(g))+abs(C(z)*C(g))+B(1))
        @test upper(defect)<=upper(polybudget)
        @test av>=0 && real(-im*Gl[e,m,1,1])>=0 && real(im*Gg[e,m,1,1])>=0
    end
    green=GreenState(GR,Gl,Gg,A,ones(3,2),scales)
    # Stored families are literal arrays with separate allocations; all six
    # nodes (including five nonzero REAL retarded nodes) enter the final map.
    fam(coef)=SelfEnergyFamily(coef.*copy(GR),coef.*copy(Gl),coef.*copy(Gg))
    total,plus,minus=fam(Float64(mid(a.c))),fam(Float64(mid(sq(a.tau)))),fam(Float64(mid(sq(a.tau))))
    checked_plain_fields((GR,Gl,Gg,A,total.Σᴿ,total.Σˡ,total.Σᵍ,
        plus.Σᴿ,plus.Σˡ,plus.Σᵍ,minus.Σᴿ,minus.Σˡ,minus.Σᵍ),
        green.condition_number,green.dyson_scale,(3,2,1,1))
    (green=green,total=total,plus=plus,minus=minus,roots=roots,stale=stale)
end

# No target Hamiltonian/residual builder on the expected side. Actual checked
# static fields only contribute bounded INPUT errors to literal h.
function hamiltonian_budget(problem,a,U,m)
    phi2=(B(1)/B(4),B(1)/B(2),B(1)/B(4))
    h=a.H0+foldl(+, (phi2[j]*B(U[j]) for j in 1:3);init=B(0))+a.t0*B((m-1)^2)
    basis=problem.basis
    projected=B(0)
    deviation=B(0)
    for j in 1:3
        actual=sq(abs(C(basis.Φ[j,1])))
        projected=projected+actual*abs(B(U[j]))
        deviation=deviation+abs(actual-phi2[j])*abs(B(U[j]))
    end
    # Independent SI expression for actual coefficient; gamma32 also covers
    # C0 eV/nm conversion and the source's short kinetic multiplication chain.
    coeff=sq(B(hbarSI))/(B(2)*B(m0SI)*B(problem.scales.E₀_J)*sq(B(problem.scales.L₀_m)))
    kinetic=coeff*sq(B(problem.grids.κ[m]))*B(real(basis.M⁻¹[1,1]))
    eta=abs(C(basis.H₀[1,1])-C(a.H0))+deviation+abs(kinetic-a.t0*B((m-1)^2))+
        abs(B(imag(basis.M⁻¹[1,1])))*coeff*sq(B(problem.grids.κ[m]))+
        gamma(32)*(abs(C(basis.H₀[1,1]))+projected+abs(kinetic))
    h,eta
end

function predicted_family(problem,state,which)
    # Exact coefficient from independently checked static T on actual input
    # side. Complex multiply implementation error stays separate (gamma16).
    tp=sq(abs(C(problem.basis.T₊[1,1])))
    tm=sq(abs(C(problem.basis.T₋[1,1])))
    q=which===:total ? tp+tm : which===:plus ? tp : tm
    X=(state.green.Gᴿ,state.green.Gˡ,state.green.Gᵍ)
    values=[[C(q)*C(x) for x in field] for field in X]
    errors=[[gamma(16)*q*abs(C(x)) for x in field] for field in X]
    values,errors
end

function norm_complex_boxes(xs)
    sqrt(foldl(+, (sq(abs(x)) for x in xs);init=B(0)))
end
function norm_positive_boxes(xs)
    sqrt(foldl(+, (sq(x) for x in xs);init=B(0)))
end

function metric_enclosures(problem,a,state,U)
    g=state.green
    expected,errors=predicted_family(problem,state,:total)
    dyson_stored=Ball[]
    dyson_candidate=Ball[]
    for m in 1:2,e in 1:3
        idx=e+3(m-1)
        h,etaH=hamiltonian_budget(problem,a,U,m)
        for (boxes,candidate) in ((dyson_stored,false),(dyson_candidate,true))
            sigma=candidate ? expected[1][idx] : C(state.total.Σᴿ[e,m,1,1])
            etaS=candidate ? errors[1][idx] : B(0)
            D0=C(B(e-2)-h)-sigma
            gg=C(g.Gᴿ[e,m,1,1])
            n0=abs(D0*gg-C(1)); d0=abs(D0)*abs(gg)+B(1)
            bn=(etaH+etaS)*abs(gg)+gamma(64)*d0
            bd=(etaH+etaS)*abs(gg)+gamma(64)*d0
            push!(boxes,ratio(n0,bn,d0,bd))
        end
    end
    # Actual independent short products, all six entries, full-array norms.
    # Candidate uncertainty includes literal TGT†'s two complex products;
    # stored uncertainty is zero because its Float64 arrays are exact inputs.
    function contraction_metric(output,sg,sgerr)
        defects=CBall[]; absolute_errors=Ball[]
        for i in eachindex(output)
            gr=C(g.Gᴿ[i]); rhs=gr*sg[i]*conj(gr)
            push!(defects,C(output[i])-rhs)
            push!(absolute_errors,sq(abs(gr))*sgerr[i]+gamma(16)*abs(rhs))
        end
        numerator=norm_complex_boxes(defects)
        normalizer=normbox(output)+B(1e-14)
        bn=norm_positive_boxes(absolute_errors)+gamma(32)*(normbox(output)+norm_complex_boxes([C(output[i])-defects[i] for i in eachindex(output)]))
        bd=gamma(32)*normalizer
        ratio(numerator,bn,normalizer,bd)
    end
    storedGamma=[C(im)*(C(state.total.Σᵍ[i])-C(state.total.Σˡ[i])) for i in eachindex(g.A)]
    candidateGamma=[C(im)*(expected[3][i]-expected[2][i]) for i in eachindex(g.A)]
    gammaError=[errors[3][i]+errors[2][i]+gamma(16)*abs(candidateGamma[i]) for i in eachindex(g.A)]
    rA=contraction_metric(g.A,storedGamma,[gamma(16)*abs(x) for x in storedGamma])
    rAc=contraction_metric(g.A,candidateGamma,gammaError)
    rK=contraction_metric(g.Gˡ,expected[2],errors[2])
    sigma_bounds=Ball[]
    for (tag,fam) in ((:total,state.total),(:plus,state.plus),(:minus,state.minus))
        pred,err=predicted_family(problem,state,tag)
        for (j,name) in enumerate((:Σᴿ,:Σˡ,:Σᵍ))
            old=getfield(fam,name)
            delta=[pred[j][i]-C(old[i]) for i in eachindex(old)]
            n0=norm_complex_boxes(delta)
            d0=norm_complex_boxes(pred[j])+B(1e-14)
            bn=norm_positive_boxes(err[j])+gamma(32)*(norm_complex_boxes(pred[j])+normbox(old))
            bd=norm_positive_boxes(err[j])+gamma(32)*d0
            push!(sigma_bounds,ratio(n0,bn,d0,bd))
        end
    end
    # Weighted raw charge, positive single-support term: no cancellation,
    # normalization, fitted lambda or hidden missing channels.
    measure=B(2)/(B(2)*a.pi)*B(problem.grids.wᴱ[2])*B(problem.grids.wᵏ[1])
    gr=C(g.Gᴿ[2,1,1,1])
    raw=measure*abs(C(im)*C(state.total.Σˡ[2,1,1,1])*(-C(1))) * sq(abs(gr))
    cand=measure*abs(expected[2][2])*sq(abs(gr))
    normed=measure*B(real(-im*g.Gˡ[2,1,1,1]))
    # Input target density conversion chain: physical sheet SI * actual L0².
    target=B(ustrip(u"m^-2",problem.physical.N_dop²ᴰ))*sq(B(problem.scales.L₀_m))
    targeterr=gamma(16)*abs(target)
    @assert lower(target)>upper(targeterr)>0
    @assert upper(abs(target-a.D))<=upper(gamma(16)*a.D)
    numbererr(x)=gamma(32)*abs(x)+measure*sq(abs(gr))*errors[2][2]
    rawbox=expand(raw,numbererr(raw)); cb=expand(cand,numbererr(cand))
    nb=expand(normed,gamma(32)*abs(normed)); tb=expand(target,targeterr)
    lambda=rawbox/tb; lc=cb/tb
    rsum=maxbox([abs(B(real(g.A[2,1,1,1]))/(B(2)*a.pi)-B(1)),B(1)])
    rsum=expand(rsum,gamma(32)*(rsum+B(1)))
    # One support node. Each real collision product has <=3 multiplications;
    # no cancellation-based relative current assertion near zero.
    A=B(real(g.A[2,1,1,1])); occ=B(real(-im*g.Gˡ[2,1,1,1]))
    fluxbound=B(2)*gamma(3)*sq(a.tau)*occ*(A-occ)/(B(16)*sq(a.pi))
    # Actual T differs from tau only at the accepted static-builder error;
    # use larger actual coefficient if needed, before any current measurement.
    qplus=sq(abs(C(problem.basis.T₊[1,1])))
    fluxbound=fluxbound*(B(1)+abs(qplus-sq(a.tau))/sq(a.tau))
    rJbound=B(2)*fluxbound/B(1e-12)
    @assert upper(rJbound)<problem_tolerances().r_J
    (rD=maxbox(dyson_stored),rDc=maxbox(dyson_candidate),rA=rA,rAc=rAc,rK=rK,
     rSigma=maxbox(sigma_bounds),raw=rawbox,candidate=cb,normalized=nb,target=tb,
     lambda=lambda,lc=lc,rlambda=abs(lambda-B(1)),rlc=abs(lc-B(1)),
     rsum=rsum,flux=fluxbound,rJ=rJbound)
end
problem_tolerances()=SolverOptions().tolerances

function poisson_enclosures(problem,a,U)
    # Independent floating stencil qL*M from source_unit_literals, not from
    # actual grid/Poisson outputs. Strict new unit chain replaces falsifiedD1.
    @assert problem.grids.wˣ==fill(a.units.h,3)
    @assert all(==(problem.profiles.Nᴰ[1]),problem.profiles.Nᴰ)
    @assert all(isfinite,problem.profiles.Nᴰ) && minimum(problem.profiles.Nᴰ)>0
    la=B(problem.scales.λ_P); nd=B(problem.profiles.Nᴰ[1])
    @assert lower(la)>0 && lower(nd)>0
    @assert upper(abs(la-a.la))<=upper(gamma(16)*a.la)
    @assert upper(abs(nd-a.D/(B(3)*a.h)))<=upper(gamma(16)*a.D/(B(3)*a.h))
    q=B(a.units.qL); v=B(a.units.v)
    # Published rounded n_f is exact data; S is its exact high-precision sum,
    # never replaced by ideal D/h. Uniform donor error affects zeta, not P*n.
    nf=B.(a.n)
    S=foldl(+,nf;init=B(0))
    meanN=S/B(3)
    b=[la*(nf[j]-nd) for j in 1:3]
    exactU=[la/(B(3)*q)*(nf[j]-meanN) for j in 1:3]
    exactZetaSqrt=la*(S-B(3)*nd)/(B(4)*q*sqrt(B(3)))
    exactZetaV=foldl(+,b;init=B(0))/(B(12)*q*v)
    neutrality=abs(S-B(3)*nd)/(B(3)*nd)
    # Exact independent Q4 top-left=M/4, v64 source inverse-sqrt coefficient.
    # Its eigenvalues are3/4 twice and +/-sqrt(3)*v; hence inverse2norm4/3.
    @assert lower(sqrt(B(3))*v)>3/4
    invQnorm=B(4)/B(3)
    Q4=[i<=3 && j<=3 ? B(a.L[i,j])/B(4) :
        (i==4 && j==4 ? B(0) : v) for i in 1:4,j in 1:4]
    rhs4=[b[1]/(B(4)*q),b[2]/(B(4)*q),b[3]/(B(4)*q),B(0)]
    zshort=gamma(32)*(la*a.D/(q*a.h)+abs(exactZetaV))
    reciprocal_v_difference=abs(exactZetaSqrt-exactZetaV)
    # Comparison uses exact v64 reference, so no ideal-v error is silently
    # inserted. The difference remains explicit and bounded independently.
    @assert upper(reciprocal_v_difference)<=upper(zshort)
    # All exact refs/caps frozen BEFORE sole tiny public solve in this helper.
    Up,zeta,rP,rneutral=solve_periodic_poisson(problem,copy(a.n))
    @test length(Up)==3 && all(isfinite,Up) && isfinite(zeta)
    rows=[q*foldl(+, (B(a.L[j,k])*B(Up[k]) for k in 1:3);init=B(0)) for j in 1:3]
    residual=[rows[j]-b[j] for j in 1:3]
    bnorm=norm_positive_boxes([abs(x) for x in b])
    rownorm=norm_positive_boxes([abs(x) for x in rows])
    # Source b formation<=3 operations; row/dot/norm<=32. qL formation was
    # separately guarded exact data, never approximated by1/h² here.
    berror=gamma(3)*bnorm
    rowerror=gamma(32)*(rownorm+bnorm)
    rounderror=rowerror+berror
    meanU=foldl(+, (B(x) for x in Up);init=B(0))/B(3)
    meanR=foldl(+,residual;init=B(0))/B(3)
    projected=[x-meanR for x in residual]
    # Literal zero-mean eigenvalue3qL and explicit gauge bound; no guessed
    # forward-solve coefficient. Formation/row error retains its h/q units.
    forward=norm_positive_boxes([abs(x) for x in projected])/(B(3)*q)+
        abs(meanU)+rounderror/(B(3)*q)
    for j in 1:3
        @test upper(abs(B(Up[j])-exactU[j]))<=upper(forward)
    end
    @test upper(abs(meanU))<=upper(forward)
    d=bnorm+la*B(1e-14)
    pbox=ratio(norm_positive_boxes([abs(x) for x in residual]),rounderror,
        d,gamma(32)*d+berror)
    check_enclosed(rP,pbox)
    @test rP<=problem_tolerances().r_P
    @assert upper(pbox)<problem_tolerances().r_P
    # Full four-row exact residual of the returned floats, without a second
    # solve. It includes source RHS formation and solver error. It must FIT
    # the retained gamma32 cap; measured error cannot enlarge acceptance.
    answer4=B.([Up...,zeta])
    R4=[foldl(+, (Q4[i,j]*answer4[j] for j in 1:4);init=B(0))-rhs4[i] for i in 1:4]
    zforward=invQnorm*norm_positive_boxes([abs(x) for x in R4])
    @assert upper(zforward)<=lower(zshort)
    @test upper(abs(B(zeta)-exactZetaV))<=upper(zforward)
    zb=expand(abs(exactZetaV),zforward) # smaller than retained regression cap
    # Actual neutrality uses weighted h dot products. The cancellation budget
    # uses ABSOLUTE terms over positive donor sheet, never relative error in0.
    donorSheet=B(3)*a.h*nd
    nnum=abs(a.h*(S-B(3)*nd))
    nerror=gamma(32)*(donorSheet+a.h*foldl(+, (abs(x) for x in nf);init=B(0)))
    nb=ratio(nnum,nerror,donorSheet,gamma(32)*donorSheet)
    @test upper(abs(neutrality-nnum/donorSheet))<=upper(gamma(32)*(B(1)+neutrality))
    check_enclosed(abs(zeta),zb)
    check_enclosed(rneutral,nb)
    @assert upper(zb)<problem_tolerances().r_ζ
    @assert upper(nb)<problem_tolerances().r_neutral
    ub=maxbox([abs(exactU[j]-B(U[j])) for j in 1:3])
    ub=expand(ub,forward+gamma(3)*(maxbox([abs(x) for x in exactU])+maxbox([abs(B(x)) for x in U])))
    (rP=pbox,rU=ub,neutral=nb,zeta=zb,Up=copy(Up),forward=forward)
end

function independent_checks(problem,a,state,bounds)
    g=state.green
    @test fieldnames(GreenState)==(:Gᴿ,:Gˡ,:Gᵍ,:A,:condition_number,:dyson_scale)
    @test fieldnames(SelfEnergyFamily)==(:Σᴿ,:Σˡ,:Σᵍ)
    @test fieldnames(FixedPointAuditCandidate)==(:scattering,:embedding,:plus,:minus,:roundoff)
    @test fieldnames(SCBAResult)==(:green,:scattering,:embedding,:embedding_plus,:embedding_minus,
        :history,:converged,:status,:quality,:restart_contract,:mixer_state)
    @test fieldnames(NEGFSolution)==(:problem,:options,:Uᴴ,:n,:scba,:outer_history,:observables,
        :report,:converged,:status)
    @test fieldnames(OuterIteration)==(:μ,:r_P,:r_U,:r_n,:r_neutral,:r_J,:r_Jchange,
        :r_population,:ζ,:J,:populations,:density)
    @test fieldnames(SCBAIteration)==(:ν,:r_D,:r_A,:r_K,:r_Σ,:r_λ,:λ,:r_PSD,:r_caus,
        :r_roundoff,:r_Jchange,:r_population,:J,:raw_charge,:target_charge,:normalized_charge,
        :lambda_change,:witness,:physical_markers)
    for e in 1:3,m in 1:2
        i=(e,m,1,1)
        gr,gl,gg,av=(X[i...] for X in (g.Gᴿ,g.Gˡ,g.Gᵍ,g.A))
        @test abs(gr-conj(gr)-(gg-gl))<=Float64(upper(gamma(16)*(B(abs(gr))+B(abs(gg))+B(abs(gl)))))
        @test av==im*(gr-conj(gr))
        @test gl+conj(gl)==gg+conj(gg)==0
        @test real(av)>=0 && real(-im*gl)>=0 && real(im*gg)>=0
        @test real(-im*gl)<=real(av)
        @test g.condition_number[e,m]==1
        for fam in (state.total,state.plus,state.minus)
            sr,sl,sg=(getfield(fam,f)[i...] for f in (:Σᴿ,:Σˡ,:Σᵍ))
            identity_bound=gamma(16)*(B(abs(sr))+B(abs(sl))+B(abs(sg)))
            @test abs(sr-conj(sr)-(sg-sl))<=Float64(upper(identity_bound))
            @test sl+conj(sl)==sg+conj(sg)==0
            @test imag(sr)<=0 && real(-im*sl)>=0 && real(im*sg)>=0
        end
        for f in (:Σᴿ,:Σˡ,:Σᵍ)
            lhs=getfield(state.total,f)[i...]
            rhs=getfield(state.plus,f)[i...]+getfield(state.minus,f)[i...]
            @test abs(lhs-rhs)<=Float64(upper(gamma(3)*(B(abs(lhs))+B(abs(rhs)))))
        end
    end
    # Supplemental metadata checks use the actual Hamiltonian only AFTER
    # literal h and etaH are frozen; target output never defines the oracle.
    inputU=state.stale ? [0.0,0.0,0.0] : copy(a.U)
    actualh=project_hamiltonians(problem,inputU)
    for m in 1:2
        h,etaH=hamiltonian_budget(problem,a,inputU,m)
        @test upper(abs(C(actualh[m,1,1])-C(h)))<=upper(etaH)
        for e in 1:3
            d=ComplexF64((e-2)-actualh[m,1,1]-state.total.Σᴿ[e,m,1,1])
            @test cond(reshape([d],1,1))==g.condition_number[e,m]
            idealSelf=C(a.c)*C(g.Gᴿ[e,m,1,1])
            deltaS=abs(C(state.total.Σᴿ[e,m,1,1])-idealSelf)
            scaleErr=etaH+deltaS+gamma(32)*(abs(C(d))+B(g.dyson_scale[e,m]))
            @test upper(abs(B(opnorm(reshape([d],1,1),Inf))-B(g.dyson_scale[e,m])))<=upper(scaleErr)
        end
    end
    # Physical independent charge/density/populations include correct E/k/spin
    # measures. The only support node gives sheet N=D and volume n=D phi²/h.
    density=ustrip.(u"m^-3",electron_density(problem,g.Gˡ)).*a.units.L0^3
    sheet=ustrip.(u"m^-2",sheet_density_matrix(problem,g.Gˡ)).*a.units.L0^2
    pop=state_populations(problem,g.Gˡ)
    number=number_functional(g.Gˡ,problem.grids,problem.physical.g_s)
    check_enclosed(number,bounds.normalized)
    chiNorm=B(0)
    chiNormInputError=B(0)
    spatialCharge=foldl(+, (a.h*B(x) for x in a.n);init=B(0))
    spatialChargeError=foldl(+, (a.h*abs(B(a.n[j])-a.nideal[j]) for j in 1:3);init=B(0))
    @test upper(abs(spatialCharge-a.D))<=upper(spatialChargeError+gamma(32)*a.D)
    for j in 1:3
        chi=sq(abs(C(problem.basis.χ[j,1])))
        ideal=a.nideal[j] # D phi²/h, independent continuous volume density
        inputError=abs(chi*a.D-ideal)
        # Includes six-entry integration and SI divide/convert/multiply chain;
        # RN(sqrt(h64))==1 need not equal exact real sqrt(h64).
        err=inputError+gamma(32)*(chi*a.D+ideal)
        check_enclosed(density[j],expand(ideal,err))
        roundedError=err+abs(B(a.n[j])-ideal)
        check_enclosed(density[j],expand(B(a.n[j]),roundedError))
        chiNorm=chiNorm+a.h*chi
        chiNormInputError=chiNormInputError+a.h*abs(chi-ideal/a.D)
    end
    @test upper(abs(chiNorm-B(1)))<=upper(chiNormInputError+gamma(32)*(chiNorm+B(1)))
    check_enclosed(real(sheet[1,1]),expand(a.D,gamma(32)*a.D))
    check_enclosed(pop.p[1],expand(B(1),gamma(32)*B(1)+abs(bounds.target-a.D)/a.D))
    @test size(sheet)==(1,1) && length(pop.p)==1
    @test bounds.rA.hi<problem_tolerances().r_A
    @test bounds.rAc.hi<problem_tolerances().r_K
    @test bounds.rK.hi<problem_tolerances().r_K
    @test bounds.rSigma.hi<problem_tolerances().r_Σ
    @test bounds.rlambda.hi<problem_tolerances().r_λ
    @test bounds.rlc.hi<problem_tolerances().r_λ
    nothing
end

function manufactured_solution(problem,a,state,stored_bounds,publishedU,poisson)
    green=deepcopy(state.green)
    total,plus,minus=deepcopy(state.total),deepcopy(state.plus),deepcopy(state.minus)
    # Dependent carrier observation from same fields, repeated identically.
    # Never an unrelated J=0; the independent zero-current expectation is the
    # absolute cancellation bound frozen in stored_bounds.
    flux=boundary_flux(problem,plus,green)
    J=eSI*Float64(ustrip(u"m^-2*s^-1",flux))
    Jbar=B(J)/B(J0ref)
    @test upper(abs(Jbar))<=upper(stored_bounds.flux*(B(1)+gamma(32)))
    populations=copy(state_populations(problem,green.Gˡ).p)
    history=SCBAIteration[]
    for nu in 1:4
        # Manufactured constant exact-algebraic sequence. Residual entries are
        # independent arithmetic upper bounds, not solver-run measurements.
        push!(history,SCBAIteration(nu,Float64(upper(stored_bounds.rD)),
            Float64(upper(stored_bounds.rA)),Float64(upper(stored_bounds.rK)),
            Float64(upper(stored_bounds.rSigma)),Float64(upper(stored_bounds.rlambda)),
            1.0,0.0,0.0,0.0,nu==1 ? Inf : 0.0,nu==1 ? Inf : 0.0,J,
            1/64,1/64,1/64,nu==1 ? NaN : 0.0,nothing,nothing))
    end
    scba=SCBAResult(green,Dict{Symbol,SelfEnergyFamily}(),total,plus,minus,history,
        true,:converged,:strictly_converged,nothing,QCLNEGF.SCBAMixerState())
    outer=OuterIteration[]
    # Identical density samples give zero change. Final r_n uses the source
    # input-period normalizer D_actual/sqrt(a.units.Lp), never sqrt(3*h).
    for mu in 1:3
        push!(outer,OuterIteration(mu,Float64(upper(poisson.rP)),
            state.stale ? maximum(abs,a.U) : 0.0,
            mu==1 ? Inf : 0.0,Float64(upper(poisson.neutral)),0.0,
            mu==1 ? Inf : 0.0,mu==1 ? Inf : 0.0,0.0,J,
            copy(populations),copy(a.n)))
    end
    options=SolverOptions() # all default scientific gates remain unchanged
    @test length(history)-1>=options.convergence.required_consecutive_scba_passes
    @test length(outer)-1>=options.convergence.required_consecutive_poisson_passes
    @test isinf(history[1].r_Jchange) && isinf(history[1].r_population)
    @test isnan(history[1].lambda_change)
    @test isinf(outer[1].r_n) && isinf(outer[1].r_Jchange) && isinf(outer[1].r_population)
    @test all(row->row.J==J,history) && all(row->row.J==J,outer)
    @test last(history).r_Σ<options.tolerances.r_Σ
    @test scba.restart_contract===nothing && scba.mixer_state.method===:linear
    @test isempty(scba.mixer_state.states) && isempty(scba.mixer_state.residuals)
    observables=Dict{Symbol,Any}(
        :fixture_kind=>"complete manufactured raw-map functional alignment",
        :history_origin=>"constant algebraic sequence; no SCBA/outer solve",
        :hartree_origin=>"literal neutral M3/h² U=lambda(n-ND), zero gauge",
        :independent_target=>1/64,
        :manufactured_scope=>"bias0/channels-off/:none/rank1/legacy bulk feedback; not S1",
        :acceptance_measurement_context=>Dict{String,Any}(
            "backend"=>"educational","fresh_candidate_required"=>true,
            "hilbert_roundoff_measured"=>false,"hilbert_applicability"=>"empty scattering set"))
    NEGFSolution(problem,options,copy(publishedU),copy(a.n),scba,outer,observables,
        ConvergenceReport(false,Dict{Symbol,Float64}(),["manufactured unassessed input"]),
        false,:manufactured_final_alignment)
end

function check_candidate(problem,state,scba)
    candidate=QCLNEGF._reference_audit_candidate(problem,scba)
    @test isempty(candidate.scattering) && isempty(scba.scattering)
    @test candidate.roundoff==0 # defined maximum of EMPTY reconstruction set
    for (tag,fam) in ((:total,candidate.embedding),(:plus,candidate.plus),(:minus,candidate.minus))
        expected,errors=predicted_family(problem,state,tag)
        for (j,f) in enumerate((:Σᴿ,:Σˡ,:Σᵍ)),i in eachindex(getfield(fam,f))
            actual=getfield(fam,f)[i]
            delta=abs(C(actual)-expected[j][i])
            @test upper(delta)<=upper(errors[j][i])
        end
    end
    arrays=(scba.green.Gᴿ,scba.green.Gˡ,scba.green.Gᵍ,scba.green.A,
        (getfield(f,k) for f in (scba.embedding,scba.embedding_plus,scba.embedding_minus,
          candidate.embedding,candidate.plus,candidate.minus) for k in (:Σᴿ,:Σˡ,:Σᵍ))...)
    checked_plain_fields(arrays,scba.green.condition_number,scba.green.dyson_scale,(3,2,1,1))
    candidate
end

failed_names(report)=Set(Symbol(first(split(message,'='))) for message in report.messages)
function assert_inventory(report,expected)
    @test all(k->haskey(report.metrics,k)&&isfinite(report.metrics[k]),GATES)
    @test !report.passed
    @test length(report.messages)==length(expected)
    @test failed_names(report)==expected
end

function check_final(solution,candidate,bounds,poisson,expected)
    # Actual public complete interface, independently frozen enclosures above.
    fresh=fresh_fixed_point_metrics(solution.problem,solution.scba,solution.Uᴴ,candidate)
    report=validate_solution(solution;candidate=candidate)
    actual_limits=QCLNEGF._solution_validation_limits(
        solution.options.tolerances,report.metrics)
    @test Set(keys(actual_limits))==Set(GATES)
    explicit=QCLNEGF._validate_solution_metrics(solution,report.metrics;
        accept_approximate_scba=false)
    @test explicit.passed==report.passed && explicit.messages==report.messages
    assert_inventory(report,expected)
    assert_inventory(explicit,expected)
    for (key,value) in fresh
        @test haskey(report.metrics,key)
        @test isequal(report.metrics[key],value)
    end
    for (key,box) in ((:r_D,bounds.rD),(:r_D_candidate,bounds.rDc),(:r_A,bounds.rA),
        (:r_A_candidate,bounds.rAc),(:r_K,bounds.rK),(:r_Σ,bounds.rSigma),
        (:raw_charge,bounds.raw),(:candidate_raw_charge,bounds.candidate),
        (:normalized_charge,bounds.normalized),(:target_charge,bounds.target),
        (:lambda,bounds.lambda),(:candidate_lambda,bounds.lc),(:r_λ,bounds.rlambda),
        (:r_lambda_candidate,bounds.rlc))
        check_enclosed(fresh[key],box)
    end
    for (key,box) in ((:r_sum,bounds.rsum),(:r_P,poisson.rP),(:r_U,poisson.rU),
        (:r_neutral,poisson.neutral),(:r_ζ,poisson.zeta))
        check_enclosed(report.metrics[key],box)
    end
    @test 0<=report.metrics[:r_J]<=upper(bounds.rJ)
    @test fresh[:r_Σ]<solution.options.tolerances.r_Σ
    @test last(solution.scba.history).r_Σ<solution.options.tolerances.r_Σ
    report
end

@testset "Final state U G alignment" begin
    # Empty metrics query only the unconditional limit inventory; this is
    # never a synthetic metric dictionary supplied for acceptance.
    base_expected=setdiff(Set(GATES),Set(CANDIDATE_GATES))
    base_actual=QCLNEGF._solution_validation_limits(
        SolverOptions().tolerances,Dict{Symbol,Float64}())
    @test length(base_expected)==38
    @test Set(keys(base_actual))==base_expected
    setprecision(BigFloat,384) do
        a=independent_literals()
        problem=manufactured_problem(3,1,:none,a)
        analytic_basis_phase(problem,a)
        @assert validate_problem(problem).passed
        @assert problem.grids.wˣ==fill(a.units.h,3)
        @assert all(isfinite,problem.profiles.Nᴰ) && minimum(problem.profiles.Nᴰ)>0
        for nd in problem.profiles.Nᴰ
            @assert upper(abs(B(nd)-a.D/(B(3)*a.h)))<=upper(gamma(16)*a.D/(B(3)*a.h))
        end
        donor=foldl(+, (B(problem.grids.wˣ[j])*B(problem.profiles.Nᴰ[j]) for j in 1:3);init=B(0))
        @assert upper(abs(donor-a.D))<=upper(gamma(32)*a.D)
        @test problem.physical.g_s==2 && problem.physical.f_ion==1
        @test fieldnames(typeof(problem))==(:physical,:numerical,:scattering,:scales,
            :grids,:profiles,:basis,:kernels,:W₊ᴱᵖ,:W₋ᴱᵖ,:W₊ᴸᴼ,:W₋ᴸᴼ,
            :energy_shift_discretization,:models)
        @test ustrip(u"K",problem.physical.Tᴸ)==ustrip(u"K",problem.physical.Tᴸᴼ)==70
        @test problem.models.electron_electron.mode===:none
        aligned=literal_state(problem,a,false)
        stale=literal_state(problem,a,true)
        # Freeze all expected metric enclosures BEFORE target fresh/final APIs.
        healthy_bounds=metric_enclosures(problem,a,aligned,a.U)
        wrongU=2 .* a.U
        wrong_bounds=metric_enclosures(problem,a,aligned,wrongU)
        stale_bounds=metric_enclosures(problem,a,stale,a.U)
        stale_stored_bounds=metric_enclosures(problem,a,stale,[0.0,0.0,0.0])
        healthy_poisson=poisson_enclosures(problem,a,a.U)
        wrong_poisson=poisson_enclosures(problem,a,wrongU)
        stale_poisson=poisson_enclosures(problem,a,a.U)
        t=SolverOptions().tolerances
        @assert upper(healthy_bounds.rD)<t.r_D && upper(healthy_bounds.rDc)<t.r_K
        @assert upper(stale_stored_bounds.rD)<t.r_D
        separation=(a.u/sqrt(a.c))/(B(2)+a.u/sqrt(a.c))
        @test lower(separation)>t.r_K
        @test lower(wrong_bounds.rD)>lower(separation)-upper(gamma(64))
        @test lower(stale_bounds.rD)>lower(separation)-upper(gamma(64))
        @test lower(wrong_bounds.rDc)>t.r_K && lower(stale_bounds.rDc)>t.r_K
        @test lower(wrong_poisson.rU)>t.r_U
        @test upper(healthy_poisson.rU)<t.r_U && upper(stale_poisson.rU)<t.r_U
        @test lower(healthy_bounds.rsum)>t.r_sum && lower(stale_bounds.rsum)>t.r_sum
        independent_checks(problem,a,aligned,healthy_bounds)
        independent_checks(problem,a,stale,stale_stored_bounds)
        healthy=manufactured_solution(problem,a,aligned,healthy_bounds,a.U,healthy_poisson)
        wrong=manufactured_solution(problem,a,aligned,healthy_bounds,wrongU,healthy_poisson)
        old=manufactured_solution(problem,a,stale,stale_stored_bounds,a.U,stale_poisson)
        hc=check_candidate(problem,aligned,healthy.scba)
        wc=check_candidate(problem,aligned,wrong.scba)
        oc=check_candidate(problem,stale,old.scba)
        # Variant copies cannot alias and silently repair a deliberately bad input.
        allarrays=Any[]
        for solution in (healthy,wrong,old)
            append!(allarrays,(solution.Uᴴ,solution.n,solution.scba.green.Gᴿ,
                solution.scba.green.Gˡ,solution.scba.green.Gᵍ,solution.scba.green.A,
                solution.scba.green.condition_number,solution.scba.green.dyson_scale))
            for fam in (solution.scba.embedding,solution.scba.embedding_plus,solution.scba.embedding_minus)
                append!(allarrays,(fam.Σᴿ,fam.Σˡ,fam.Σᵍ))
            end
            for row in solution.outer_history
                append!(allarrays,(row.populations,row.density))
            end
        end
        for candidate in (hc,wc,oc),fam in (candidate.embedding,candidate.plus,candidate.minus)
            append!(allarrays,(fam.Σᴿ,fam.Σˡ,fam.Σᵍ))
        end
        @test all(!Base.mightalias(allarrays[i],allarrays[j]) for i in eachindex(allarrays) for j in (i+1):length(allarrays))
        hr=check_final(healthy,hc,healthy_bounds,healthy_poisson,Set((:r_sum,)))
        wr=check_final(wrong,wc,wrong_bounds,wrong_poisson,Set((:r_D,:r_D_candidate,:r_U,:r_sum)))
        sr=check_final(old,oc,stale_bounds,stale_poisson,Set((:r_D,:r_D_candidate,:r_sum)))
        # Negative copies are controls of the test contract, never input to a
        # production evaluator: old/history metrics cannot satisfy the frozen
        # bad-state oracle; predicted U cannot satisfy retained-input rU.
        for (bad,bnd) in ((wr,wrong_bounds),(sr,stale_bounds))
            substituted=copy(bad.metrics)
            substituted[:r_D]=hr.metrics[:r_D]
            substituted[:r_D_candidate]=hr.metrics[:r_D_candidate]
            @test !contains(bnd.rD,substituted[:r_D])
            @test !contains(bnd.rDc,substituted[:r_D_candidate])
        end
        repaired=copy(wrong.Uᴴ)
        repaired.=healthy_poisson.Up
        @test repaired!=wrong.Uᴴ && wrong.Uᴴ==wrongU
        @test !contains(wrong_poisson.rU,hr.metrics[:r_U])
        @test healthy.Uᴴ==old.Uᴴ==a.U
        @test healthy.scba.green.Gᴿ!=old.scba.green.Gᴿ
        @test !healthy.converged && !wrong.converged && !old.converged
        @test all(s->s.status===:manufactured_final_alignment,(healthy,wrong,old))
        # A functional aligned state still fails spectral quadrature. It is
        # not a verified discretization, solver convergence or experimental result.
    end
end
end
