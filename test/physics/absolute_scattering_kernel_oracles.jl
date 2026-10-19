module AbsoluteScatteringKernelOracles
include("../support/common.jl")

# DEV02 O08-05 / DEV07, accepted D1 (2026-10-19).
# A manufactured Nb=1 coefficient fixture, not the S1 scientific calculation.
# Independent SI covariance derivations and applicability: see the DEV07
# subsection of docs/src/development/review_oracle_contracts.md.
# No expected value calls CODATA, _scaled_physics, _K₀, kernels or reducers.

# Small outward interval implementation: EACH primitive is rounded outward.
# Operation inventory, source inspection of kernels.jl at Core 7324d98.
# Counts are conservative REAL operation ceilings, not dynamic profiling.
# Unary conjugation/sign changes/copies/indexing are exact; complex multiply6,
# complex addition2, real*complex and complex/real2. Reused expressions may
# be charged again, so this counts work rather than relying on reuse.
# SI constants/parameter conversion, K0 preparation and conditioned band-jump
# subtraction: <=128 equivalent relative-error factors for every channel.
# SI allowance expansion: <=16 literal-to-Float64 roundings, <=32 unit
# conversion/scale factors (including K0 four operations), <=32 material
# energy construction/contrast factors (subtraction condition9 charged),
# <=32 prefactor operations, <=16 spare counted conversion factors =128.
# This concerns the fixed canonical units above, not arbitrary Unitful paths.
# This includes the subtraction condition9 for 2.5tau-2tau; normal positive
# physical inputs are checked separately. Grid errors use g and are propagated
# separately, so eigensolver/grid builder operations are NOT silently counted
# as a 1024-step eigensolver: independent residual/gap certificates replace it.
#
# Independent CP-CODE conservative ceilings (E1 accounting erratum):
# LO: SI128 +9 envelope/phase terms*24
#     +4 angles*(q12 +3 qz blocks*24 +scale12) +accumulation16 =744.
# A qz block counts three squares and two adds (denominator5), four
# factor operations and complex products/accumulation10: 19<24.
# IFR: SI128 +3 interfaces*32
#      +4 angles*(q12 +spectrum16 +3 vertex contractions*12 +scale12)
#      +accumulation16 =544.
# Acoustic: SI128 +3 quartic terms*32 +prefactor16 +scale12
#           +extra accumulation16 =268.
# Alloy: SI128 +3 variance/quartic terms*36 +prefactor16 +scale12
#        +extra accumulation16 =280.
# These independently checked ceilings are all <992<1024. exp/sin/cos
# errors have their explicit MPFR certificate; sqrt assumes IEEE correct
# rounding. All products are normal and finite for the fixed SI literals.
# Evaluating a whole nonmonotone expression twice in global rounding modes
# would not enclose it. All use is inside setprecision(BigFloat,256).
struct IV
    lo::BigFloat
    hi::BigFloat
end
down(f) = setrounding(f, BigFloat, RoundDown)
up(f) = setrounding(f, BigFloat, RoundUp)
IV(x::AbstractString) = IV(down(() -> parse(BigFloat,x)), up(() -> parse(BigFloat,x)))
IV(x::Real) = IV(down(() -> BigFloat(x)), up(() -> BigFloat(x)))
Base.:+(a::IV,b::IV) = IV(down(() -> a.lo+b.lo), up(() -> a.hi+b.hi))
Base.:-(a::IV,b::IV) = IV(down(() -> a.lo-b.hi), up(() -> a.hi-b.lo))
Base.:-(a::IV) = IV(-a.hi,-a.lo)
function Base.:*(a::IV,b::IV)
    lows = [down(() -> x*y) for x in (a.lo,a.hi), y in (b.lo,b.hi)]
    highs = [up(() -> x*y) for x in (a.lo,a.hi), y in (b.lo,b.hi)]
    IV(minimum(lows),maximum(highs))
end
function Base.inv(a::IV)
    (a.lo>0 || a.hi<0) || error("interval division through zero")
    IV(down(() -> inv(a.hi)),up(() -> inv(a.lo)))
end
Base.:/(a::IV,b::IV) = a*inv(b)
for op in (:+,:-,:*,:/)
    @eval Base.$op(a::IV,b::Real) = Base.$op(a,IV(b))
    @eval Base.$op(a::Real,b::IV) = Base.$op(IV(a),b)
end
function sq(a::IV)
    lo = a.lo<=0<=a.hi ? BigFloat(0) : min(down(() -> a.lo*a.lo),down(() -> a.hi*a.hi))
    IV(lo,max(up(() -> a.lo*a.lo),up(() -> a.hi*a.hi)))
end
function Base.sqrt(a::IV)
    a.lo>=0 || error("negative interval square root")
    IV(down(() -> sqrt(a.lo)),up(() -> sqrt(a.hi)))
end
Base.exp(a::IV) = IV(down(() -> exp(a.lo)),up(() -> exp(a.hi)))
pi_interval() = IV(down(() -> BigFloat(pi)),up(() -> BigFloat(pi)))
absupper(a::IV) = max(abs(a.lo),abs(a.hi))
mid(a::IV) = (a.lo+a.hi)/2
sumiv(xs) = foldl(+,xs;init=IV(0))

function arithmetic_bounds()
    u=IV("1.1102230246251565404236316680908203125e-16")
    gamma(n)=n*u/(1-n*u)
    g,a=gamma(32),gamma(1024)
    beta=g/(1+sqrt(1-g))+sqrt(IV(2))*g*sqrt(IV(29)/2)/((IV(3)/2)*sqrt(1-g))
    h=2beta+sq(beta)
    eQ=g+(1+g)*(2h+3sq(h))/(IV(13)/21)
    bI=beta+6g/(1-4g)+g*(1+beta)
    hI=2bI+sq(bI)
    eI=(2hI*(IV(207)/112)+sq(hI)*(IV(13)/2))/(IV(2329)/3136)
    dq=128g
    dx=dq/4+16g
    dIFR=(1+a)*(1+eI)*exp(dx)*(1+2u)-1
    eta=(5pi_interval()/2)*(2g+sq(g))+2sqrt(IV(2))*u
    eW=3h*(1+g)+g
    ferr=eW+(1+eW)*eta+a*(1+eW)*(1+eta)
    eF2=(2ferr+sq(ferr))/(IV(169)/441)
    dD=dq+(sq(pi_interval())+1)*(2g+sq(g))+a
    @test dD.hi<1
    dLO=(1+a)*(1+g)*(1+eF2)/(1-dD)-1
    delta=Dict(:LO=>dLO,:IFR=>dIFR,:acoustic=>(1+a)*(1+eQ)-1,
               :alloy=>(1+a)*(1+eQ)-1)
    @test all(d -> 0<d.lo<=d.hi<IV(1//8).lo,values(delta))
    (;g,beta,ferr,delta)
end

function literals()
    e,hbar,kB,m0,eps0,ell,mr=IV.(
        ("1.602176634e-19","1.054571817e-34","1.380649e-23",
         "9.1093837139e-31","8.8541878128e-12","1e-8","0.067"))
    tau=sq(hbar)/(2m0*mr*sq(ell))
    (;e,hbar,kB,m0,eps0,ell,mr,tau)
end

function carrier_inputs(;T=70.0,qs=1e8,height=1e-10,x=0.25,disabled=false,
                        missing_energy=false,missing_volume=false)
    c=literals()
    tau_eV=Float64(mid(c.tau/c.e))
    layers=[Layer(d=1e-8u"m",material=:manufactured,x_Al=x,
        Eᶜ=(v*tau_eV)*u"eV",mᶻᵣ=0.067,m_parallelᵣ=0.067,εᵣ=16.0,doped=true)
        for v in (2.0,2.5,0.5)]
    # PhysicalParameters has a checked positional constructor. Canonical
    # units here match its concrete field types; no unchecked carrier structs.
    p=PhysicalParameters(layers,3e16u"m^-2",1.0,0.0u"V/m",0.0u"m",0.0u"eV",
        T*u"K",70.0u"K",16.0,8.0,0.05u"eV",1e8u"m^-1",qs*u"m^-1",
        height*u"m",1e-8u"m",4.0u"eV",4000.0u"kg/m^3",4000.0u"m/s",
        (disabled || missing_energy) ? nothing : 0.6u"eV",
        (disabled || missing_volume) ? nothing : 3.125e-29u"m^3",2,
        [0.0u"m",1e-8u"m",2e-8u"m"])
    n=NumericalParameters(N_z=3,N_b=1,P_basis=1,E_min=-0.1u"eV",E_max=0.1u"eV",
        N_E=3,M_E=0.05u"eV",k_max=2e8u"m^-1",N_k=2,N_φ=4,
        qz_max=(Float64(pi)/1e-8)*u"m^-1",N_qz=3,η_seed=0.001u"eV")
    p,n
end

# Every Float64 input/unit conversion is checked against decimal/rational SI
# literals before it can enter a kernel. This is separate from reference RHS.
function scalar_certificate(actual,ref::IV,g::IV)
    @test isfinite(actual)
    error=absupper(IV(actual)-ref)
    @test error <= (g*IV(max(BigFloat(1),absupper(ref)))).lo
end
function vector_certificate(actual,refs,g)
    @test length(actual)==length(refs)
    for (x,r) in zip(actual,refs)
        scalar_certificate(x,r,g)
    end
end
function si_certificate(actual,ref::IV,g::IV)
    @test isfinite(actual)
    if iszero(ref.lo) && iszero(ref.hi)
        @test iszero(actual)
    else
        @test absupper(IV(actual)/ref-1)<=g.lo
    end
end
function input_certificate(p,n,s,E0,L0,b)
    c=literals()
    si_certificate(ustrip(u"J",uconvert(u"J",p.layers[1].Eᶜ)),2c.tau,b.g)
    for (layer,v) in zip(p.layers,(IV(2),IV(5)/2,IV(1)/2))
        si_certificate(ustrip(u"J",uconvert(u"J",layer.Eᶜ)),v*c.tau,b.g)
        si_certificate(ustrip(u"m",layer.d),c.ell,b.g)
        si_certificate(layer.mᶻᵣ,c.mr,b.g)
        si_certificate(layer.m_parallelᵣ,c.mr,b.g)
        @test layer.εᵣ==16 && layer.doped
    end
    for (value,unit,ref) in (
        (p.N_dop²ᴰ,u"m^-2",IV("3e16")),(p.F_bias,u"V/m",IV(0)),
        (p.z₀,u"m",IV(0)),(p.E_ref,u"eV",IV(0)),
        (p.Tᴸ,u"K",IV(70)),(p.Tᴸᴼ,u"K",IV(70)),
        (p.ħωᴸᴼ,u"J",IV(".05")*c.e),(p.q_s,u"m^-1",IV("1e8")),
        (p.qᴸᴼ_s,u"m^-1",IV("1e8")),
        (p.Δᴵᶠᴿ,u"m",IV("1e-10")),
        (p.Λᴵᶠᴿ,u"m",c.ell),(p.Ξ,u"J",4c.e),
        (p.ρ_m,u"kg/m^3",IV(4000)),(p.v_s,u"m/s",IV(4000)))
        si_certificate(ustrip(unit,uconvert(unit,value)),ref,b.g)
    end
    @test p.f_ion==1 && p.ε_s==16 && p.ε_∞==8 && p.g_s==2
    if p.ΔV_alloy!==nothing
        si_certificate(ustrip(u"J",uconvert(u"J",p.ΔV_alloy)),IV(".6")*c.e,b.g)
        si_certificate(ustrip(u"m^3",p.Ω₀),IV("3.125e-29"),b.g)
    end
    for (z,r) in zip(p.interfaces,[IV(0),c.ell,2c.ell])
        si_certificate(ustrip(u"m",z),r,b.g)
    end
    si_certificate(s.E₀_eV,E0,b.g)
    si_certificate(s.E₀_J,E0*c.e,b.g)
    si_certificate(s.L₀_m,L0,b.g)
    si_certificate(ustrip(u"m^-1",n.qz_max),pi_interval()/c.ell,b.g)
    for (value,unit,ref) in ((n.E_min,u"eV",-IV(".1")),
        (n.E_max,u"eV",IV(".1")),(n.M_E,u"eV",IV(".05")),
        (n.k_max,u"m^-1",IV("2e8")),(n.η_seed,u"eV",IV(".001")))
        si_certificate(ustrip(unit,value),ref,b.g)
    end
    @test (n.N_z,n.N_b,n.P_basis,n.N_E,n.N_k,n.N_φ,n.N_qz)==(3,1,1,3,2,4,3)
end

# Grid certificates use dimensionless coordinates scaled to cell length;
# exact zeros do not get silently replaced by missing data or an atol.
function geometry_certificate(g,pf,basis,lambda,E0,b)
    c=literals(); piI=pi_interval()
    vector_certificate(g.x,[lambda*(IV(j)-IV(1)/2) for j=1:3],b.g)
    vector_certificate(g.wˣ,fill(lambda,3),b.g)
    vector_certificate(g.κ,[IV(0),2/lambda],b.g)
    vector_certificate(g.wᵏ,[1/(4piI*sq(lambda)),3/(4piI*sq(lambda))],b.g)
    vector_certificate(g.qᶻ,[-piI/lambda,IV(0),piI/lambda],b.g)
    vector_certificate(g.wᑫᶻ,[piI/(2lambda),piI/lambda,piI/(2lambda)],b.g)
    vector_certificate(g.φ,[IV(0),piI/2,piI,3piI/2],b.g)
    vector_certificate(g.ε,[-IV(".1")/E0,IV(0),IV(".1")/E0],b.g)
    vector_certificate(g.wᴱ,[IV(".05")/E0,IV(".1")/E0,IV(".05")/E0],b.g)
    @test iszero(g.κ[1]) && iszero(g.qᶻ[2])
    @test g.trusted_energy==[false,true,false]
    vector_certificate(pf.Eᶜ,[v*c.tau/(E0*c.e) for v in (IV(2),IV(5)/2,IV(1)/2)],b.g)
    vector_certificate(pf.mᶻᵣ,fill(c.mr,3),b.g)
    vector_certificate(pf.m_parallelᵣ,fill(c.mr,3),b.g)
    vector_certificate(pf.εᵣ,fill(IV(16),3),b.g)
    vector_certificate(pf.Nᴰ,fill(inv(lambda*sq(lambda)),3),b.g)
    @test pf.layer_index==[1,2,3]
    @test pf.x_Al==fill(pf.x_Al[1],3)
    @test basis.localization===:none && size(basis.χ)==(3,1)
    # Independently evaluate ||(A-2I)c|| and normalization in interval arithmetic.
    # Complex components are kept separate, including a possible common phase.
    cr=[sqrt(lambda)*IV(real(z)) for z in basis.χ[:,1]]
    ci=[sqrt(lambda)*IV(imag(z)) for z in basis.χ[:,1]]
    Aminus=[2 -1 0; -1 5//2 -1; 0 -1 1//2]
    norm2=sumiv([sq(cr[j])+sq(ci[j]) for j=1:3])
    residual2=sumiv([sq(sumiv([Aminus[i,j]*v[j] for j=1:3]))
                     for i=1:3 for v in (cr,ci)])
    @test absupper(norm2-1)<=b.g.lo
    @test sqrt(residual2).hi<=(b.g*sqrt(IV(29)/2)).lo
    # Explicit alignment with c3/|c3|; c3 is bounded away from zero by gap.
    phase_norm=sqrt(sq(cr[3])+sq(ci[3]))
    @test phase_norm.lo>0
    pr,pii=cr[3]/phase_norm,ci[3]/phase_norm
    v=[IV(j)/sqrt(IV(21)) for j in (1,2,4)]
    error2=sumiv([sq(cr[j]-v[j]*pr)+sq(ci[j]-v[j]*pii) for j=1:3])
    @test sqrt(error2).hi<=b.beta.lo
end

# MPFR correctly rounded scalar functions at exact represented Float64 inputs.
# The test bounds error by ONE adjacent Float64 spacing, not arbitrary eps.
function math_certificate(f,x::Float64,actual::Float64=f(x))
    @test isfinite(x) && isfinite(actual)
    X=BigFloat(x) # exact binary input (precision=256)
    lo=down(() -> f(X)); hi=up(() -> f(X))
    nearest_lo=Float64(lo); nearest_hi=Float64(hi)
    @test nearest_lo==nearest_hi # certify the correctly rounded Float64 value
    r=nearest_lo
    ulp=min((IV(nextfloat(r))-IV(r)).lo,(IV(r)-IV(prevfloat(r))).lo)
    @test absupper(IV(actual)-IV(lo,hi))<=ulp
end
function library_certificates(g,L0)
    for phi in g.φ
        math_certificate(cos,phi)
    end
    for qz in g.qᶻ,x in g.x
        # exp(im*qz*x) may use complex exp rather than scalar sin/cos entrypoints.
        theta=qz*x
        z=exp(im*qz*x)
        math_certificate(cos,theta,real(z))
        math_certificate(sin,theta,imag(z))
    end
    for k in g.κ,kp in g.κ,phi in g.φ
        qbar=sqrt(max(k^2+kp^2-2k*kp*cos(phi),0.0))
        q=qbar/L0
        arg=-q^2*1e-8^2/4
        math_certificate(exp,arg)
    end
    for arg in (0.0,-1.0,-2.0,-4.0)
        math_certificate(exp,arg)
    end
end

function expected(channel,E0,L0;screen=1,T=70)
    c=literals(); piI=pi_interval()
    K0=sq(E0*c.e)*sq(L0)
    function f(v)
        if channel===:LO
            I=c.ell/2*(1/(IV(v)+screen^2)+(IV(169)/441)/(IV(v)+screen^2+sq(piI)))
            return sq(c.e)*(IV(".05")*c.e)/(32c.eps0)*I/K0
        elseif channel===:IFR
            return piI*sq(IV("1e-10"))*sq(c.tau)*(IV(2329)/3136)*exp(-IV(v)/4)/K0
        elseif channel===:acoustic
            return sq(4c.e)*c.kB*T/(4000*IV(4000)*4000)*(IV(13)/(21c.ell))/K0
        elseif channel===:alloy
            return IV("3.125e-29")*sq(IV(".6")*c.e)*(IV(13)/(112c.ell))/K0
        end
        error("unknown channel")
    end
    [f(0) f(4); f(4) (f(0)+2f(8)+f(16))/4]
end
function accepts(z,R::IV,delta::IV)
    isfinite(z) || return false
    R.lo>0 || error("positive oracle required")
    lower=R*(1-delta); upper=R*(1+delta)
    BigFloat(real(z))>=lower.lo && BigFloat(real(z))<=upper.hi &&
        BigFloat(abs(imag(z)))<=(R*delta).hi
end
function tensor_assertions(K,R,delta)
    @test size(K)==(2,2,1,1,1,1)
    for i=1:2,j=1:2
        @test accepts(K[i,j,1,1,1,1],R[i,j],delta)
        @test !accepts(2K[i,j,1,1,1,1],R[i,j],delta)
    end
    @test K./maximum(abs,K)==(2K)./maximum(abs,2K)
end

@testset "Absolute raw scattering covariance amplitudes (Nb=1)" begin
    setprecision(BigFloat,256) do
        b=arithmetic_bounds()
        @test_throws ArgumentError carrier_inputs(missing_energy=true)
        @test_throws ArgumentError carrier_inputs(missing_volume=true)
        raw_cases=Dict{Symbol,Any}[]
        for (E0,L0) in ((IV(".1"),IV("1e-8")),(IV(".2"),IV("2e-8")))
            p,n=carrier_inputs()
            s=ScaleSystem(E₀=Float64(mid(E0))*u"eV",L₀=Float64(mid(L0))*u"m")
            g=build_grids(p,n,s); pf=build_profiles(p,n,s,g)
            basis=build_localized_basis(p,n,s,g,pf;localization=:none)
            input_certificate(p,n,s,E0,L0,b)
            geometry_certificate(g,pf,basis,literals().ell/L0,E0,b)
            @test pf.x_Al==fill(0.25,3)
            library_certificates(g,s.L₀_m)
            LO,F=lo_phonon_kernel(p,s,g,basis)
            @test size(F)==(3,1,1)
            for (z,w) in zip(vec(F),(complex(0,-13//21),complex(1,0),complex(0,13//21)))
                @test absupper(IV(real(z))-IV(real(w)))<=b.ferr.lo
                @test absupper(IV(imag(z))-IV(imag(w)))<=b.ferr.lo
            end
            raw=Dict(:LO=>LO,:IFR=>interface_roughness_kernel(p,s,g,basis),
                     :acoustic=>acoustic_kernel(p,s,g,basis),:alloy=>alloy_kernel(p,s,g,pf,basis))
            for channel in (:LO,:IFR,:acoustic,:alloy)
                R=expected(channel,E0,L0)
                tensor_assertions(raw[channel],R,b.delta[channel])
                # Named microscopic wrong-operator factors: primary normalization
                # mistakes must fail independently of normalized shape.
                factors=channel===:LO ? (2Float64(pi),8.8541878128e-12) :
                        channel===:IFR ? (inv(Float64(pi)),2313/2329,2345/2329) :
                        channel===:acoustic ? (0.5,1e-8) : (4.0,0.5,4/3)
                for factor in factors
                    @test !accepts(factor*raw[channel][1,2,1,1,1,1],R[1,2],b.delta[channel])
                end
            end
            # Gaussian exponent -v/2, nearest-cell interface amplitudes:
            @test !accepts(raw[:IFR][1,2,1,1,1,1]*exp(-1.0),expected(:IFR,E0,L0)[1,2],b.delta[:IFR])
            nearest=(IV(1)/14)*(IV(1)/14)+(IV(2)/21)*(IV(2)/21)+(IV(32)/21)*(IV(32)/21)
            @test !accepts(raw[:IFR][1,2,1,1,1,1]*Float64(mid(nearest/(IV(2329)/3136))),expected(:IFR,E0,L0)[1,2],b.delta[:IFR])
            p2,_=carrier_inputs(qs=2e8)
            LO2,_=lo_phonon_kernel(p2,s,g,basis)
            tensor_assertions(LO2,expected(:LO,E0,L0;screen=2),b.delta[:LO])
            p35,_=carrier_inputs(T=35.0)
            tensor_assertions(acoustic_kernel(p35,s,g,basis),expected(:acoustic,E0,L0;T=35),b.delta[:acoustic])
            pzero,_=carrier_inputs(height=0.0)
            @test all(iszero,interface_roughness_kernel(pzero,s,g,basis))
            for x in (0.0,1.0)
                px,_=carrier_inputs(x=x)
                pfx=build_profiles(px,n,s,g)
                @test all(iszero,alloy_kernel(px,s,g,pfx,basis))
            end
            pd,_=carrier_inputs(disabled=true)
            @test all(iszero,alloy_kernel(pd,s,g,pf,basis))
            push!(raw_cases,raw)
        end
        # Supplementary units control AFTER the separate absolute comparisons.
        for channel in (:LO,:IFR,:acoustic,:alloy),i=1:2,j=1:2
            A=raw_cases[1][channel][i,j,1,1,1,1]
            B=raw_cases[2][channel][i,j,1,1,1,1]
            d=b.delta[channel]
            @test abs(BigFloat(real(A))-16BigFloat(real(B))) <=
                  (expected(channel,IV(".1"),IV("1e-8"))[i,j]*2d).hi
        end
    end
end
end # module
