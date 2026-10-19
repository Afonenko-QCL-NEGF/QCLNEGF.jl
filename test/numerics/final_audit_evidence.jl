module FinalAuditEvidence
include("../support/common.jl")
const BN = QCLNEGF
const QN = QCLNEGF.QCLNumerics
const QR = QCLNEGF.QCLReferenceOperators
# Manufactured algebra only; this file does not certify shared map kernels or
# a device solution. Expected values never call a product residual/number/solve.
const u = 2.0^-53
const δ = Float64(1e-14)
gamma(n) = n*u/(1-n*u)
const Bmatrix = 2.0^24*gamma(4096)
const COMPONENTS = ("R", "L", "G")
scalar(x) = reshape(ComplexF64[x], 1,1,1,1)
family(r=0,l=0,g=0) = SelfEnergyFamily(scalar(r),scalar(l),scalar(g))
state(r=0,l=0,a=0) = GreenState(scalar(r),scalar(l),scalar(l-im*a),scalar(a),ones(1,1),ones(1,1))
value(m) = m["value"]
measured(x) = Dict{String,Any}("status"=>"measured","value"=>Float64(x),"nonfinite"=>nothing,"reason"=>nothing)

@testset "Same-pass independent Dyson witnesses and nonfinite inventory" begin
    gr = reshape(ComplexF64[1,1],2,1,1,1)
    sigma = reshape(ComplexF64[1,-10],2,1,1,1) # E=h=0 => D=(-1,10)
    h = zeros(ComplexF64,1,1,1)
    observed = QN._dyson_residual_evidence([0.0,0.0],h,sigma,gr)
    @test observed.ratio == 1
    @test abs(observed.global_abs_F-sqrt(85.0)) <= 2u*sqrt(85.0)
    @test observed.max_abs_F == 9
    rw, aw = observed.ratio_witness, observed.absolute_witness
    @test (rw.e_index,rw.k_index) == (1,1)
    @test (aw.e_index,aw.k_index) == (2,1)
    @test (rw.numerator_F,rw.D_F,rw.GR_F,rw.identity_F,rw.denominator) == (2,1,1,1,2)
    @test abs(aw.ratio-9/11) <= gamma(2)*abs(9/11)
    @test aw.denominator == 11
    @test observed.global_abs_F != rw.numerator_F # witness-only mutation
    @test (aw.e_index,aw.k_index) != (rw.e_index,rw.k_index)
    @test BN._dyson_residual([0.0,0.0],h,sigma,gr) === observed.ratio
    healthy = QN._dyson_residual_evidence([1.0,1.0],h,zeros(ComplexF64,size(sigma)),gr)
    @test healthy.ratio == healthy.global_abs_F == healthy.max_abs_F == 0
    @test (healthy.ratio_witness.e_index,healthy.ratio_witness.k_index) == (1,1)
    @test (healthy.absolute_witness.e_index,healthy.absolute_witness.k_index) == (1,1)
    for bad in (NaN,Inf,-Inf)
        corrupt = copy(sigma); corrupt[2,1,1,1] = bad
        result = QN._dyson_residual_evidence([0.0,0.0],h,corrupt,gr)
        @test !isfinite(result.ratio)
        @test (result.ratio_witness.e_index,result.ratio_witness.k_index) == (2,1)
        @test (result.absolute_witness.e_index,result.absolute_witness.k_index) == (2,1)
        @test QN._audit_witness_record(result.ratio_witness)["status"] == "error"
    end
end

@testset "Literal near-zero additive floor and measured-zero controls" begin
    a = 2.0^-60
    # BigFloat represents the actual Float64 floor, not decimal 1e-14.
    expected = setprecision(BigFloat,256) do
        Float64(BigFloat(a)/(BigFloat(a)+BigFloat(δ)))
    end
    bound = gamma(4)*abs(expected)
    A = QN._spectral_residual_evidence(state(0,0,a),family())
    K = QN._keldysh_residual_evidence(state(0,im*a,0),family())
    S = QN._selfenergy_residual_evidence(family(),family(a)).rows[1]
    for raw in (A,K,S)
        @test raw.numerator_F == raw.scale_F == a
        @test raw.denominator == a+δ
        @test abs(raw.ratio-expected) <= bound
        @test abs(a/max(a,δ)-expected) > bound
        @test abs(1.0-expected) > bound # omitted floor mutant
    end
    for raw in (QN._spectral_residual_evidence(state(),family()),
                QN._keldysh_residual_evidence(state(),family()),
                QN._selfenergy_residual_evidence(family(),family()).rows[1])
        @test raw.numerator_F == raw.scale_F == raw.ratio == 0
        @test raw.denominator == δ
        @test QN._audit_norm_record(raw)["ratio"] == measured(0)
    end
end

# Literal determinant/adjugate and 2x2 scalar sums. No inv, solve, norm,
# number_functional or tested reducer is used for these expected quantities.
mm(a,b) = [a[i,1]*b[1,j]+a[i,2]*b[2,j] for i=1:2,j=1:2]
adj(a) = [conj(a[j,i]) for i=1:2,j=1:2]
frob(a) = sqrt(sum(real(z)^2+imag(z)^2 for z in a))
function inverse2(a)
    d = a[1,1]*a[2,2]-a[1,2]*a[2,1]
    return [a[2,2]/d -a[1,2]/d; -a[2,1]/d a[1,1]/d], d
end

function matrix_fixture()
    energy, we, wk = [-1/2,3/4],[1/2,3/2],[1/4,3/4]
    hs = ([0 (1+im)/8; (1-im)/8 1/4],
          [1/8 (1+2im)/16; (1-2im)/16 3/8])
    Γ = [1/2 im/8; -im/8 3/4]
    f = [1/4 3/8; 1/2 5/8]
    gr,gl,gg,aa,sr,sl,sg = (zeros(ComplexF64,2,2,2,2) for _=1:7)
    h = zeros(ComplexF64,2,2,2)
    dets = ComplexF64[]
    raw_number, spectral2, spectral_scale2, occupied_scale2 = 0.0,0.0,0.0,0.0
    for e=1:2,k=1:2
        h[k,:,:] .= hs[k]
        d = [energy[e]*(i==j)-hs[k][i,j]+im*Γ[i,j]/2 for i=1:2,j=1:2]
        r, determinant = inverse2(d)
        push!(dets,determinant)
        l = im*f[e,k].*Γ; g = -im*(1-f[e,k]).*Γ
        lesser = mm(mm(r,l),adj(r))
        spectral = im.*(r-adj(r))
        gr[e,k,:,:] .= r; gl[e,k,:,:] .= lesser
        gg[e,k,:,:] .= mm(mm(r,g),adj(r)); aa[e,k,:,:] .= spectral
        sr[e,k,:,:] .= -im.*Γ/2; sl[e,k,:,:] .= l; sg[e,k,:,:] .= g
        raw_number += 2/(2π)*we[e]*wk[k]*real(-im*(lesser[1,1]+lesser[2,2]))
        spectral2 += frob(spectral-mm(mm(r,Γ),adj(r)))^2
        spectral_scale2 += frob(spectral)^2
        occupied_scale2 += frob(lesser)^2
    end
    return (; energy,we,wk,h,dets,raw_number,
            spectral_numerator=sqrt(spectral2),spectral_scale=sqrt(spectral_scale2),
            occupied_scale=sqrt(occupied_scale2),
            green=GreenState(gr,gl,gg,aa,ones(2,2),ones(2,2)),
            sigma=SelfEnergyFamily(sr,sl,sg))
end

@testset "Frozen complex matrix, weighted number and unnormalized candidate" begin
    x = matrix_fixture()
    @test x.dets == [65/256-23im/64,7/16-7im/16,65/256+27im/64,1/8+11im/32]
    @test Bmatrix < 7.63e-6
    @test x.spectral_scale > 1/2 && x.occupied_scale > 1/8 && x.raw_number > 1/16
    D = QN._dyson_residual_evidence(x.energy,x.h,x.sigma.Σᴿ,x.green.Gᴿ)
    A = QN._spectral_residual_evidence(x.green,x.sigma)
    K = QN._keldysh_residual_evidence(x.green,x.sigma)
    @test D.global_abs_F <= Bmatrix
    @test abs(A.numerator_F-x.spectral_numerator) <= Bmatrix
    @test abs(A.scale_F-x.spectral_scale) <= Bmatrix
    @test K.numerator_F <= Bmatrix
    @test abs(K.scale_F-x.occupied_scale) <= Bmatrix
    grids = ModelGrids([0.0],[1.0],x.energy,x.we,trues(2),[0.0,1.0],x.wk,[0.0],[1.0],[0.0])
    observed_number = number_functional(x.green.Gˡ,grids,2)
    @test abs(observed_number-x.raw_number) <= Bmatrix
    q = 9/8
    candidate = SelfEnergyFamily(copy(x.sigma.Σᴿ),q.*x.sigma.Σˡ,
                                q.*x.sigma.Σˡ.+x.sigma.Σᵍ.-x.sigma.Σˡ)
    fresh = QN._keldysh_residual_evidence(x.green,candidate)
    expected_num = (q-1)*x.occupied_scale
    expected_ratio = expected_num/(x.occupied_scale+δ)
    @test abs(fresh.numerator_F-expected_num) <= Bmatrix
    @test abs(fresh.ratio-expected_ratio) <= Bmatrix
    candidate_raw = number_functional(keldysh_green(x.green.Gᴿ,candidate.Σˡ),grids,2)
    @test abs(candidate_raw-q*x.raw_number) <= Bmatrix
    @test abs(candidate_raw/x.raw_number-q) <= Bmatrix
    # Actual observed mutants: dividing the candidate by q masks the defect;
    # coherence/ordinary-transpose/wrong-weight changes violate literal oracles.
    divided = SelfEnergyFamily(copy(candidate.Σᴿ),candidate.Σˡ./q,copy(candidate.Σᵍ))
    masked = QN._keldysh_residual_evidence(x.green,divided)
    @test abs(masked.ratio-expected_ratio) > Bmatrix
    mutated_gl = copy(x.green.Gˡ); mutated_gl[1,1,1,2] += 1
    mutated = GreenState(x.green.Gᴿ,mutated_gl,x.green.Gᵍ,x.green.A,ones(2,2),ones(2,2))
    @test QN._keldysh_residual_evidence(mutated,x.sigma).numerator_F > Bmatrix
    wrong_transpose = zeros(ComplexF64,2,2,2,2)
    wrong_number = 0.0
    for e=1:2,k=1:2
        r = x.green.Gᴿ[e,k,:,:]; l = x.sigma.Σˡ[e,k,:,:]
        wrong = mm(mm(r,l),permutedims(r))
        wrong_transpose[e,k,:,:] .= wrong
        wrong_number += 2/(2π)*x.we[e]*x.wk[k]*real(-im*(wrong[1,1]+wrong[2,2]))
    end
    @test abs(number_functional(wrong_transpose,grids,2)-x.raw_number) > Bmatrix
    @test abs(wrong_number-x.raw_number) > Bmatrix
    wrong_grids = ModelGrids([0.0],[1.0],x.energy,reverse(x.we),trues(2),[0.0,1.0],x.wk,[0.0],[1.0],[0.0])
    @test abs(number_functional(x.green.Gˡ,wrong_grids,2)-x.raw_number) > Bmatrix
    @test abs(number_functional(x.green.Gˡ,grids,1)-x.raw_number) > Bmatrix
    sigma = QN._selfenergy_residual_evidence(x.sigma,candidate)
    literal_difference = frob(vec((q-1).*x.sigma.Σˡ))
    literal_scale = frob(vec(q.*x.sigma.Σˡ))
    @test abs(sigma.rows[2].numerator_F-literal_difference) <= Bmatrix
    @test abs(sigma.rows[2].ratio-literal_difference/(literal_scale+δ)) <= Bmatrix
    mixed = SelfEnergyFamily(copy(x.sigma.Σᴿ),x.sigma.Σˡ.+2.0^-20 .*(candidate.Σˡ.-x.sigma.Σˡ),copy(candidate.Σᵍ))
    @test abs(QN._selfenergy_residual_evidence(x.sigma,mixed).rows[2].numerator_F-literal_difference) > Bmatrix
end

function poisson_fixture()
    h, er, donor, a = 2.0,4.0,1.0,1/16
    grids = ModelGrids([1.0,3.0,5.0],fill(h,3),[0.0],[1.0],trues(1),[0.0],[1.0],[0.0],[1.0],[0.0])
    profiles = MaterialProfiles(zeros(3),ones(3),ones(3),fill(er,3),fill(donor,3),zeros(3),ones(Int,3))
    scales = ScaleSystem()
    n = donor .+ a.*[1.0,-2.0,1.0]
    expected = scales.λ_P*a*h^2/(3er).*[1.0,-2.0,1.0]
    return (; grids,profiles,scales,n,expected)
end

# Directed scalar interval graph retained from accepted O08-03 D2 §5.
# Expected L/b/Q4 never call the product builder, norm, ratio or solver.
struct PBall
    lo::BigFloat
    hi::BigFloat
end
pdown(f)=setrounding(f,BigFloat,RoundDown)
pup(f)=setrounding(f,BigFloat,RoundUp)
PB(x::Real)=PBall(BigFloat(x),BigFloat(x))
Base.:+(a::PBall,b::PBall)=PBall(pdown(()->a.lo+b.lo),pup(()->a.hi+b.hi))
Base.:-(a::PBall,b::PBall)=PBall(pdown(()->a.lo-b.hi),pup(()->a.hi-b.lo))
Base.:*(a::PBall,b::PBall)=PBall(minimum(pdown(()->x*y) for x in (a.lo,a.hi),y in (b.lo,b.hi)),
                                      maximum(pup(()->x*y) for x in (a.lo,a.hi),y in (b.lo,b.hi)))
function Base.:/(a::PBall,b::PBall)
    @assert b.lo>0 || b.hi<0
    a*PBall(pdown(()->inv(b.hi)),pup(()->inv(b.lo)))
end
Base.sqrt(a::PBall)=(@assert a.lo>=0; PBall(pdown(()->sqrt(a.lo)),pup(()->sqrt(a.hi))))
function Base.abs(a::PBall)
    a.lo>=0 && return a
    a.hi<=0 && return PBall(-a.hi,-a.lo)
    PBall(BigFloat(0),max(-a.lo,a.hi))
end
pcontains(a::PBall,x::Real)=a.lo<=BigFloat(x)<=a.hi
pgamma(n)=PB(n)*PB(u)/(PB(1)-PB(n)*PB(u))
psum(xs)=foldl(+,xs;init=PB(0))
pnorm(xs)=sqrt(psum(abs(x)*abs(x) for x in xs))
pmax(xs)=PBall(maximum(x.lo for x in xs),maximum(x.hi for x in xs))
pexpand(x,r)=PBall(max(BigFloat(0),pdown(()->x.lo-r.hi)),pup(()->x.hi+r.hi))
function pratio(n,bn,d,bd)
    @assert d.lo>bd.hi>=0
    PBall(max(BigFloat(0),pdown(()->(n.lo-bn.hi)/(d.hi+bd.hi))),
          pup(()->(n.hi+bn.hi)/(d.lo-bd.hi)))
end

function literal_poisson_reference(profiles,grids,scales,n)
    h,er,nd = grids.wˣ[1],profiles.εᵣ[1],profiles.Nᴰ[1]
    @assert length(n)==3 && grids.wˣ==fill(h,3)
    @assert profiles.εᵣ==fill(er,3) && profiles.Nᴰ==fill(nd,3)
    @assert er in (1.0,4.0) && nd>0 && h>0
    q64 = er/(h*h) # independent floating input-side coefficient, no matrix helper
    @assert all(sum(abs,row)==4q64 for row in ([2q64,-q64,-q64],[-q64,2q64,-q64],[-q64,-q64,2q64]))
    q,la,hb,donor = PB(q64),PB(scales.λ_P),PB(h),PB(nd)
    nf=PB.(n); S=psum(nf); meanN=S/PB(3)
    b=[la*(nf[j]-donor) for j=1:3]
    exactU=[la/(PB(3)*q)*(nf[j]-meanN) for j=1:3]
    v64=0x1.279a74590331dp-1 # accepted RN(1/RN(sqrt(3)))
    v=PB(v64)
    @assert (sqrt(PB(3))*v).lo>3/4
    M=[2 -1 -1;-1 2 -1;-1 -1 2]
    Q4=[i<=3 && j<=3 ? PB(M[i,j])/PB(4) :
        (i==4 && j==4 ? PB(0) : v) for i=1:4,j=1:4]
    rhs4=[b[1]/(PB(4)*q),b[2]/(PB(4)*q),b[3]/(PB(4)*q),PB(0)]
    exactZ=psum(b)/(PB(12)*q*v)
    donorSheet=PB(3)*hb*donor
    # Accepted Bz_short=γ32*(la*D/(q*h)+abs(exactZ)); D=actual donor sheet.
    short=pgamma(32)*(la*donorSheet/(q*hb)+abs(exactZ))
    return (;q,la,hb,donor,nf,S,b,exactU,exactZ,M,Q4,rhs4,short,donorSheet)
end

function literal_poisson_enclosure(ref,result,publishedU)
    Up,zeta,_,_=result.solution
    q,la,b=ref.q,ref.la,ref.b
    rows=[q*psum(PB(ref.M[j,k])*PB(Up[k]) for k=1:3) for j=1:3]
    residual=[rows[j]-b[j] for j=1:3]
    bnorm=pnorm(b); rownorm=pnorm(rows)
    berror=pgamma(3)*bnorm
    rounderror=pgamma(32)*(rownorm+bnorm)+berror
    meanU=psum(PB.(Up))/PB(3); meanR=psum(residual)/PB(3)
    forward=pnorm([x-meanR for x in residual])/(PB(3)*q)+abs(meanU)+rounderror/(PB(3)*q)
    answer4=PB.([Up...,zeta])
    R4=[psum(ref.Q4[i,j]*answer4[j] for j=1:4)-ref.rhs4[i] for i=1:4]
    qforward=PB(4)/PB(3)*pnorm(R4) # literal inverse2norm 4/3, no second solve
    den=bnorm+la*PB(δ)
    numerator=pexpand(pnorm(residual),rounderror)
    rhs=pexpand(bnorm,pgamma(32)*bnorm+berror)
    denominator=pexpand(den,pgamma(32)*den+berror)
    pbox=pratio(pnorm(residual),rounderror,den,pgamma(32)*den+berror)
    signed=ref.hb*(PB(3)*ref.donor-ref.S)
    nerror=pgamma(32)*(ref.donorSheet+ref.hb*psum(abs.(ref.nf)))
    neutral=pratio(abs(signed),nerror,ref.donorSheet,pgamma(32)*ref.donorSheet)
    rawU=pmax([abs(ref.exactU[j]-PB(publishedU[j])) for j=1:3])
    rawU=pexpand(rawU,forward+pgamma(3)*(pmax(abs.(ref.exactU))+pmax(abs.(PB.(publishedU)))))
    return (;forward,qforward,meanU,numerator,rhs,denominator,pbox,signed,nerror,neutral,rawU)
end

function literal_poisson_matches(ref,result,publishedU)
    Up,zeta,rp,rn=result.solution; d=result.diagnostics
    all(isfinite,Up) && isfinite(zeta) || return false
    box=literal_poisson_enclosure(ref,result,publishedU)
    # Measured residual cannot inflate its own permissible error: both bounds
    # must fit the retained gamma32 cap before ANY forward assertion can pass.
    box.forward.hi<=ref.short.lo && box.qforward.hi<=ref.short.lo || return false
    all(abs(PB(Up[j])-ref.exactU[j]).hi<=box.forward.hi for j=1:3) || return false
    abs(box.meanU).hi<=box.qforward.hi || return false
    abs(PB(zeta)-ref.exactZ).hi<=box.qforward.hi || return false
    all((pcontains(box.numerator,d.residual_numerator_F),pcontains(box.rhs,d.rhs_F),
         pcontains(box.denominator,d.denominator),pcontains(box.pbox,rp),
         pcontains(box.pbox,d.residual_ratio),pcontains(box.neutral,rn),
         pcontains(box.neutral,d.neutrality_ratio))) || return false
    pcontains(pexpand(ref.donorSheet,pgamma(32)*ref.donorSheet),d.target_sheet) || return false
    abs(PB(d.neutrality_signed)-box.signed).hi<=box.nerror.hi || return false
    pcontains(pexpand(abs(box.signed),box.nerror),d.neutrality_abs) || return false
    d.additive_floor==Float64(ref.la.lo)*δ && d.border_scale==4*Float64(ref.q.lo) || return false
    d.zeta_signed==zeta || return false
    return true
end

@testset "One owning Poisson solve with independent literal residual/gauge enclosure" begin
    setprecision(BigFloat,384) do
        x=poisson_fixture(); published=[1/8,-1/4,1/8]
        ref=literal_poisson_reference(x.profiles,x.grids,x.scales,x.n) # before solve
        result=QR._solve_periodic_poisson_with_diagnostics(x.profiles,x.grids,x.scales,x.n)
        @test literal_poisson_matches(ref,result,published)
        @test solve_periodic_poisson(x.profiles,x.grids,x.scales,x.n)==result.solution
        @test length(result.solution)==4
        # Actual observed-path mutants; expected ref is immutable.
        for badUp in (-result.solution[1],result.solution[1]./4,1.1 .*result.solution[1],
                      result.solution[1].+1/16)
            mutated=(;solution=(badUp,result.solution[2:end]...),diagnostics=result.diagnostics)
            @test !literal_poisson_matches(ref,mutated,published)
        end
        altered=QR._solve_periodic_poisson_with_diagnostics(x.profiles,x.grids,x.scales,x.n.+[1/16,0,-1/16])
        @test !literal_poisson_matches(ref,altered,published)
        neutralref=literal_poisson_reference(x.profiles,x.grids,x.scales,x.profiles.Nᴰ)
        neutral=QR._solve_periodic_poisson_with_diagnostics(x.profiles,x.grids,x.scales,x.profiles.Nᴰ)
        @test literal_poisson_matches(neutralref,neutral,published)
        @test neutral.diagnostics.denominator==x.scales.λ_P*δ
        @test neutral.diagnostics.rhs_F==neutral.diagnostics.residual_numerator_F==0
        # Positive raw residual makes factor-two norm/ratio tampering observable
        # beyond the gamma32 cancellation allowance of the neutral branch.
        nonneutral_n=x.n.+1/8
        nonneutralref=literal_poisson_reference(x.profiles,x.grids,x.scales,nonneutral_n)
        nonneutral=QR._solve_periodic_poisson_with_diagnostics(x.profiles,x.grids,x.scales,nonneutral_n)
        @test literal_poisson_matches(nonneutralref,nonneutral,published)
        @test nonneutral.diagnostics.neutrality_signed==-3/4
        @test nonneutral.diagnostics.neutrality_ratio==1/8
        for replacement in (
            (;residual_numerator_F=2*nonneutral.diagnostics.residual_numerator_F),
            (;residual_ratio=2*nonneutral.diagnostics.residual_ratio),
            (;residual_numerator_F=2*nonneutral.diagnostics.residual_numerator_F,
              residual_ratio=2*nonneutral.diagnostics.residual_ratio),
        )
            wrongd=merge(nonneutral.diagnostics,replacement)
            @test !literal_poisson_matches(nonneutralref,(;solution=nonneutral.solution,diagnostics=wrongd),published)
        end
    end
end

function record_fixture(; enabled=Symbol[])
    x = matrix_fixture()
    grids = ModelGrids([0.0],[1.0],x.energy,x.we,trues(2),[0.0,1.0],x.wk,[0.0],[1.0],[0.0])
    empty_family = SelfEnergyFamily(zeros(ComplexF64,2,2,2,2),zeros(ComplexF64,2,2,2,2),zeros(ComplexF64,2,2,2,2))
    physical = reference_parameters(V_period=0.0u"mV",Tᴸ=70u"K")
    options = SolverOptions()
    scba = (; green=x.green,scattering=Dict{Symbol,SelfEnergyFamily}(name=>empty_family for name in enabled),
             embedding=x.sigma,embedding_plus=empty_family,embedding_minus=empty_family,
             history=SCBAIteration[],status=:max_iterations)
    problem = (; kernels=(; enabled),models=PhysicalModelOptions(),scales=ScaleSystem(),grids,physical)
    observables = Dict{Symbol,Any}()
    solution = (; scba,problem,outer_history=OuterIteration[],observables)
    sigma = Dict{Symbol,Any}(name=>QN._selfenergy_residual_evidence(empty_family,empty_family)
                            for name in (enabled...,:embedding_total,:embedding_plus,:embedding_minus))
    evidence = (; dyson=QN._dyson_residual_evidence(x.energy,x.h,x.sigma.Σᴿ,x.green.Gᴿ),
                 spectral=QN._spectral_residual_evidence(x.green,x.sigma),
                 keldysh=QN._keldysh_residual_evidence(x.green,x.sigma),sigma)
    metrics = Dict{Symbol,Float64}(:target_charge=>x.raw_number,:raw_charge=>x.raw_number,
        :normalized_charge=>x.raw_number,:lambda=>1,:r_λ=>0,:r_charge_normalized=>0,
        :candidate_raw_charge=>x.raw_number,:candidate_lambda=>1,:r_lambda_candidate=>0,
        :lambda_change=>NaN,:r_Σ=>0,:r_U=>0,:r_D_candidate=>evidence.dyson.ratio,
        :r_A_candidate=>evidence.spectral.ratio,:r_K=>evidence.keldysh.ratio)
    p = poisson_fixture()
    poisson = QR._solve_periodic_poisson_with_diagnostics(p.profiles,p.grids,p.scales,p.n).diagnostics
    return (; solution,evidence,metrics,poisson)
end

# This literal predicate reads each required status/value/tag/reason. Mutating
# an individual M must fail it even if overall partial_error is unchanged.
function overflow_mapping_matches(record)
    dyson = record["dyson"]
    finite_error_zero = Dict("status"=>"error","value"=>0.0,
                             "nonfinite"=>nothing,"reason"=>"nonfinite_measurement")
    inf_error = Dict("status"=>"error","value"=>nothing,
                     "nonfinite"=>"+inf","reason"=>"nonfinite_measurement")
    valid_finite(m) = m["status"] == "measured" && m["value"] isa Float64 &&
                      isfinite(m["value"]) && m["nonfinite"] === m["reason"] === nothing
    record["record_status"] == "partial_error" && record["reason"] == "measurement_error" || return false
    dyson["ratio"] == finite_error_zero || return false
    all(valid_finite(dyson[n]) for n in ("global_abs_F","max_abs_F")) || return false
    for name in ("ratio_witness","absolute_witness")
        w = dyson[name]
        w["status"] == "error" && w["e_index"] == w["k_index"] == 1 || return false
        w["ratio"] == finite_error_zero || return false
        w["D_F"] == inf_error && w["denominator"] == inf_error || return false
        all(valid_finite(w[n]) for n in ("numerator_F","GR_F","identity_F")) || return false
    end
    return true
end

@testset "Approved D2 finite ratio0 with overflowed normalization" begin
    c, g = Float64(1e308), Float64(1e-308)
    gr = reshape(ComplexF64[g,0,0,g],1,1,2,2)
    sigma = fill(ComplexF64(-c),1,1,2,2)
    h = zeros(ComplexF64,1,2,2)
    @test all(isfinite,gr) && all(isfinite,sigma)
    raw = QN._dyson_residual_evidence([0.0],h,sigma,gr)
    old_ratio = BN._dyson_residual([0.0],h,sigma,gr)
    # Frozen phenotype is required. A different runtime outcome is a failure,
    # never permission to change constants, arithmetic or the criterion.
    @test raw.ratio === old_ratio === 0.0
    @test !signbit(raw.ratio)
    @test isfinite(raw.global_abs_F) && isfinite(raw.max_abs_F)
    p = Float64(c*g)
    expected_a = sqrt(2*(p-1)^2+2*p^2)
    bound = gamma(64)*max(1,expected_a)
    @test abs(raw.global_abs_F-expected_a) <= bound
    @test abs(raw.max_abs_F-expected_a) <= bound
    for w in (raw.ratio_witness,raw.absolute_witness)
        @test w.D_F === Inf && w.denominator === Inf
        @test isfinite(w.GR_F) && w.GR_F > 0
        @test w.identity_F == sqrt(2.0)
        @test isfinite(w.numerator_F) && w.numerator_F > 0
        @test abs(w.numerator_F-expected_a) <= bound
        @test w.ratio === 0.0
    end
    x = record_fixture()
    green = GreenState(gr,zeros(ComplexF64,1,1,2,2),zeros(ComplexF64,1,1,2,2),
                       zeros(ComplexF64,1,1,2,2),ones(1,1),ones(1,1))
    solution = merge(x.solution,(; scba=merge(x.solution.scba,(; green))))
    evidence = merge(x.evidence,(; dyson=raw))
    metrics = copy(x.metrics); metrics[:r_D_candidate] = old_ratio
    record = QN._final_fixed_point_record(solution,evidence,metrics,x.poisson)
    @test record["shape"] == Dict("NE"=>1,"Nk"=>1,"Nb"=>2)
    @test overflow_mapping_matches(record)
    @test record["dyson"]["ratio"]["value"] === old_ratio
    for (name,wraw) in (("ratio_witness",raw.ratio_witness),("absolute_witness",raw.absolute_witness))
        w = record["dyson"][name]
        for field in (:numerator_F,:GR_F,:identity_F)
            @test w[String(field)] == measured(getfield(wraw,field))
        end
    end
    for name in ("global_abs_F","max_abs_F")
        @test abs(value(record["dyson"][name])-expected_a) <= bound
    end
    # Actual payload mutants, inspected by exactly the healthy assertion.
    ratio_paths = (("dyson","ratio"),("dyson","ratio_witness","ratio"),
                   ("dyson","absolute_witness","ratio"))
    for path in ratio_paths
        for mutation in (:measured,:forged_tag,:missing_value)
            mutant = deepcopy(record)
            parent = mutant
            for key in path[1:end-1]; parent = parent[key]; end
            m = parent[last(path)]
            if mutation === :measured
                m["status"] = "measured"; m["reason"] = nothing
            elseif mutation === :forged_tag
                m["nonfinite"] = "+inf"
            else
                m["value"] = nothing
            end
            @test mutant["record_status"] == "partial_error"
            @test !overflow_mapping_matches(mutant)
        end
    end
    for name in ("global_abs_F","max_abs_F")
        mutant = deepcopy(record)
        mutant["dyson"][name]["status"] = "error"
        mutant["dyson"][name]["reason"] = "nonfinite_measurement"
        @test !overflow_mapping_matches(mutant)
    end
    for name in ("ratio_witness","absolute_witness")
        mutant = deepcopy(record)
        mutant["dyson"][name]["numerator_F"]["status"] = "error"
        @test !overflow_mapping_matches(mutant)
        for field in ("D_F","denominator")
            mutant = deepcopy(record)
            mutant["dyson"][name][field] = measured(0)
            @test !overflow_mapping_matches(mutant)
        end
    end
    # Same-shape healthy measured zero has a finite normalization.
    identity_gr = reshape(ComplexF64[1,0,0,1],1,1,2,2)
    healthy_raw = QN._dyson_residual_evidence([1.0],h,zeros(ComplexF64,1,1,2,2),identity_gr)
    healthy_green = GreenState(identity_gr,green.Gˡ,green.Gᵍ,green.A,ones(1,1),ones(1,1))
    healthy_solution = merge(solution,(; scba=merge(solution.scba,(; green=healthy_green))))
    healthy = QN._final_fixed_point_record(healthy_solution,merge(evidence,(; dyson=healthy_raw)),metrics,x.poisson)
    @test healthy["record_status"] == "measured"
    @test healthy["dyson"]["ratio"] == measured(0)
    @test healthy["dyson"]["global_abs_F"] == measured(0)
    @test healthy["dyson"]["max_abs_F"] == measured(0)
    for name in ("ratio_witness","absolute_witness")
        w = healthy["dyson"][name]
        @test w["status"] == "measured" && w["e_index"] == w["k_index"] == 1
        @test w["ratio"] == measured(0)
        @test w["numerator_F"] == measured(0)
        for field in ("D_F","GR_F","identity_F","denominator")
            @test w[field]["status"] == "measured"
            @test isfinite(w[field]["value"]) && w[field]["value"] > 0
            @test w[field]["nonfinite"] === w[field]["reason"] === nothing
        end
    end
end

@testset "Bounded Core template, explicit unavailable and signed/nonfinite charge" begin
    enabled = [:alloy,:LO,:IFR,:impurity,:acoustic,:electron_electron] # canonicalization tested
    x = record_fixture(; enabled)
    record = QN._final_fixed_point_record(x.solution,x.evidence,x.metrics,x.poisson)
    @test Set(keys(record)) == Set(("schema","record_status","reason","detail","detail_truncated",
        "shape","units","context","history","dyson","spectral","keldysh","selfenergy","charge","poisson"))
    @test record["schema"] == "qcl-negf-final-fixed-point-audit-v1"
    @test !haskey(record,"state_identity") # internal stage; Runner must bind before format acceptance
    @test record["record_status"] == "measured" && record["reason"] === nothing
    @test record["context"]["state_alignment"] == "unknown"
    @test record["context"]["alignment_reason"] == "producer_precondition_unverified"
    @test record["context"]["independent_physics_oracle"] === false
    @test record["context"]["map_kernels_shared_with_solver"] === nothing
    @test record["context"]["unknown_fields"] == ["backend","energy_shift","hilbert","contraction",
        "kernel_build","embedding","embedding_periods","map_kernels_shared_with_solver"]
    @test record["history"]["previous_lambda"]["reason"] == "history_unavailable"
    @test record["charge"]["lambda_change"]["status"] == "unavailable"
    expected_families = ["LO","acoustic","impurity","IFR","alloy","electron_electron",
                         "embedding_total","embedding_plus","embedding_minus"]
    @test record["context"]["enabled_scattering"] == expected_families[1:6]
    rows = record["selfenergy"]["rows"]
    @test length(rows) == 27
    @test [(row["family"],row["component"]) for row in rows] ==
          [(f,c) for f in expected_families for c in COMPONENTS]
    @test record["selfenergy"]["winner"] == Dict("family"=>"LO","component"=>"R")
    for (section,key,metric) in (("dyson","ratio",:r_D_candidate),("spectral","ratio",:r_A_candidate),
                                ("keldysh","ratio",:r_K),("selfenergy","max_ratio",:r_Σ),
                                ("charge","candidate_constraint",:r_lambda_candidate),("poisson","rU_max",:r_U))
        @test value(record[section][key]) === x.metrics[metric]
    end
    for (inner,previous_status,previous_reason) in (
        ([(; ν=7,λ=0.875)],"not_applicable","no_previous_inner_iteration"),
        ([(; ν=7,λ=0.875),(; ν=11,λ=1.0)],"measured",nothing),
    )
        contextual = merge(x.solution,(; scba=merge(x.solution.scba,(; history=inner))))
        traced = QN._final_fixed_point_record(contextual,x.evidence,x.metrics,x.poisson)
        @test traced["history"]["inner_last"] == inner[end].ν
        @test traced["history"]["previous_lambda"]["status"] == previous_status
        @test traced["history"]["previous_lambda"]["reason"] == previous_reason
        @test traced["record_status"] == "measured"
        if length(inner) == 2
            @test traced["history"]["inner_previous"] == 7
            @test traced["charge"]["lambda_change"] == measured(1/8)
        else
            @test traced["history"]["inner_previous"] === nothing
            @test traced["charge"]["lambda_change"]["status"] == "not_applicable"
        end
    end
    # Core reports only the actual producer selectors, with no default injection.
    x.solution.observables[:acceptance_measurement_context] = Dict(
        "backend"=>"production","final_state_alignment"=>"producer_aligned",
        "map_kernels_shared_with_solver"=>true)
    x.solution.observables[:restart_contract] = Dict("algorithms"=>Dict(
        "energy_shift"=>"sparse_plan","hilbert"=>"fft","contraction"=>"dense_blas",
        "kernel_build"=>"tabulated","embedding"=>"finite_chain","embedding_periods"=>4))
    executed = QN._final_fixed_point_record(x.solution,x.evidence,x.metrics,x.poisson)
    @test executed["context"]["status"] == "known"
    @test isempty(executed["context"]["unknown_fields"])
    @test executed["context"]["embedding_periods"] == 4
    @test executed["context"]["state_alignment"] == "producer_aligned"
    @test executed["context"]["map_kernels_shared_with_solver"] === true
    # Complete fresh-dependent unavailable inventory; no measured zero substitutes.
    absent = QN._final_fixed_point_record(x.solution,nothing,x.metrics,x.poisson; audit_failure="unavailable "*repeat("я",300))
    @test absent["record_status"] == "unavailable" && absent["reason"] == "map_unavailable"
    @test ncodeunits(absent["detail"]) <= 512 && absent["detail_truncated"]
    missing = Any[absent["dyson"][key] for key in ("ratio","global_abs_F","max_abs_F")]
    for key in ("ratio_witness","absolute_witness")
        witness = absent["dyson"][key]
        @test witness["status"] == "unavailable"
        @test witness["e_index"] === witness["k_index"] === nothing
        append!(missing,[witness[n] for n in ("numerator_F","D_F","GR_F","identity_F","denominator","ratio")])
    end
    for section in ("spectral","keldysh")
        append!(missing,[absent[section][n] for n in ("numerator_F","scale_F","denominator","ratio")])
    end
    @test length(absent["selfenergy"]["rows"]) == 27
    @test absent["selfenergy"]["winner"] === nothing
    push!(missing,absent["selfenergy"]["max_ratio"])
    for row in absent["selfenergy"]["rows"], n in ("numerator_F","scale_F","denominator","ratio")
        push!(missing,row[n])
    end
    append!(missing,[absent["charge"][n] for n in ("candidate_raw_sheet","candidate_lambda","candidate_constraint")])
    @test all(m == Dict("status"=>"unavailable","value"=>nothing,"nonfinite"=>nothing,"reason"=>"map_unavailable") for m in missing)
    @test absent["charge"]["stored_raw_sheet"]["status"] == "unavailable"
    @test absent["charge"]["target_sheet"]["status"] == "measured"
    @test absent["poisson"]["residual_ratio"]["status"] == "measured"
    @test_throws ArgumentError charge_constraint_diagnostics(1,0,1)
    @test_throws ArgumentError charge_constraint_diagnostics(1,-1,1)
    @test_throws ArgumentError charge_constraint_diagnostics(0,1,1)
    for target in (0.0,-1.0)
        bad = copy(x.metrics); bad[:target_charge] = target
        invalid = QN._final_fixed_point_record(x.solution,x.evidence,bad,x.poisson)
        @test invalid["record_status"] == "partial_error"
        @test invalid["charge"]["target_sheet"]["value"] === target
        @test invalid["charge"]["target_sheet"]["reason"] == "nonpositive_target"
        @test invalid["charge"]["stored_lambda"]["reason"] == "invalid_prerequisite"
    end
    bad = copy(x.metrics); bad[:candidate_raw_charge] = -x.metrics[:target_charge]
    bad[:candidate_lambda] = -1; bad[:r_lambda_candidate] = 2
    negative = QN._final_fixed_point_record(x.solution,x.evidence,bad,x.poisson)
    @test negative["charge"]["candidate_raw_sheet"]["value"] < 0
    @test negative["charge"]["candidate_lambda"] == measured(-1)
    @test negative["charge"]["candidate_constraint"] == measured(2)
    for (nonfinite,tag) in ((NaN,"nan"),(Inf,"+inf"),(-Inf,"-inf"))
        bad = copy(x.metrics); bad[:candidate_raw_charge] = nonfinite
        invalid = QN._final_fixed_point_record(x.solution,x.evidence,bad,x.poisson)
        @test invalid["record_status"] == "partial_error"
        raw = invalid["charge"]["candidate_raw_sheet"]
        @test raw["status"] == "error" && raw["value"] === nothing && raw["nonfinite"] == tag
        for key in ("candidate_lambda","candidate_constraint")
            @test invalid["charge"][key] == Dict("status"=>"error","value"=>nothing,
                                                  "nonfinite"=>nothing,"reason"=>"invalid_prerequisite")
        end
    end
end


@testset "Public validation retains failed outcome, raw rU and compatible assessment transport" begin
    setprecision(BigFloat,384) do
        ref = reference_parameters(V_period=0.0u"mV",Tᴸ=70u"K")
        layers = Layer[Layer(d=30.0u"nm",material=:manufactured_uniform,x_Al=0.0,
            Eᶜ=0.0u"eV",mᶻᵣ=0.067,m_parallelᵣ=0.067,εᵣ=1.0,doped=true)]
        physical = PhysicalParameters((name === :layers ? layers : getfield(ref,name)
                                       for name in fieldnames(PhysicalParameters))...)
        numerical = NumericalParameters(N_z=3,N_b=1,P_basis=1,N_E=3,N_k=2,N_φ=4,N_qz=3,
            E_min=-0.1u"eV",E_max=0.1u"eV",M_E=0.05u"eV",k_max=1e8u"m^-1",
            qz_max=1e8u"m^-1",η_seed=0.001u"eV")
        problem = build_problem(; physical,numerical,
            scattering=ScatteringOptions(LO=false,acoustic=false,impurity=false,IFR=false,alloy=false),
            scales=ScaleSystem(),validate_static=false)
        gr,gl,gg,aa,sr,sl,sg = (zeros(ComplexF64,3,2,1,1) for _=1:7)
        for e=1:3,k=1:2
            r = 1/(problem.grids.ε[e]+im)
            gr[e,k,1,1] = r; aa[e,k,1,1] = 2abs2(r)
            gl[e,k,1,1] = im*abs2(r)/2; gg[e,k,1,1] = -3im*abs2(r)/2
            sr[e,k,1,1] = -im; sl[e,k,1,1] = im/2; sg[e,k,1,1] = -3im/2
        end
        green = GreenState(gr,gl,gg,aa,ones(3,2),ones(3,2))
        embedding = SelfEnergyFamily(sr,sl,sg)
        empty_family = SelfEnergyFamily(zeros(ComplexF64,3,2,1,1),zeros(ComplexF64,3,2,1,1),zeros(ComplexF64,3,2,1,1))
        scba = SCBAResult(green,Dict{Symbol,SelfEnergyFamily}(),embedding,empty_family,empty_family,
                          SCBAIteration[],false,:max_iterations,:unresolved)
        observables = Dict{Symbol,Any}()
        options = SolverOptions()
        # Accepted O08-03 actual source unit chain, frozen before product calls.
        nm=1.0e-9; expected_h=((30.0*nm)/(10.0*nm))/3.0
        @test expected_h==0x1.0000000000001p+0
        @test problem.grids.wˣ==fill(expected_h,3)
        @test problem.profiles.εᵣ==ones(3)
        @test 1.0/(expected_h*expected_h)==0x1.ffffffffffffcp-1
        nd=problem.profiles.Nᴰ[1]
        @test isfinite(nd) && nd>0 && problem.profiles.Nᴰ==fill(nd,3)
        density=nd .+ (nd/16).*[1.0,-2.0,1.0]
        publishedU=[1/8,-1/4,1/8]
        pref=literal_poisson_reference(problem.profiles,problem.grids,problem.scales,density)
        # One tiny direct measurement before the independent public validation;
        # no outer/SCBA loop, expected graph uses only literals and exact input data.
        pactual=QR._solve_periodic_poisson_with_diagnostics(problem,density)
        @test literal_poisson_matches(pref,pactual,publishedU)
        pbox=literal_poisson_enclosure(pref,pactual,publishedU)
        solution = NEGFSolution(problem,options,publishedU,density,scba,OuterIteration[],
            observables,ConvergenceReport(false,Dict{Symbol,Float64}(),String[]),false,:max_iterations)
        report = validate_solution(solution; audit_failure="manufactured-unavailable")
        @test !report.passed
        @test all(isinf(report.metrics[key]) for key in (:r_K,:r_Σ,:r_λ))
        @test any(occursin("manufactured-unavailable",message) for message in report.messages)
        record = observables[:final_fixed_point_audit]
        @test record["record_status"] == "unavailable"
        @test record["context"]["state_alignment"] == "unknown"
        @test pcontains(pbox.rawU,report.metrics[:r_U])
        @test pcontains(pbox.rawU,value(record["poisson"]["rU_max"]))
        @test value(record["poisson"]["rU_max"])===report.metrics[:r_U]
        @test pcontains(pbox.numerator,value(record["poisson"]["residual_numerator_F"]))
        @test pcontains(pbox.rhs,value(record["poisson"]["rhs_F"]))
        @test pcontains(pbox.pbox,report.metrics[:r_P])
        mixed_rU=maximum(abs.(0.125 .*(pactual.solution[1].-publishedU)))
        @test !pcontains(pbox.rawU,mixed_rU)
        changed_density=density.+(nd/16).*[1.0,-2.0,1.0]
        changed = NEGFSolution(problem,options,publishedU,changed_density,scba,OuterIteration[],
            Dict{Symbol,Any}(),solution.report,false,:max_iterations)
        changed_report=validate_solution(changed; audit_failure="manufactured-unavailable")
        @test !pcontains(pbox.rawU,changed_report.metrics[:r_U])
        final = NEGFSolution(problem,options,solution.Uᴴ,solution.n,scba,OuterIteration[],
                             observables,report,false,:validation_failed)
        assessment = solution_scientific_assessment(final)
        @test assessment["final_fixed_point_audit"] === record
        @test assessment["scientific_accepted"] === false
        @test assessment["stationary_candidate_accepted"] === false
        @test all(limits[key] == options.tolerances.r_K for key in (:r_D_candidate,:r_A_candidate)
                  for limits in (BN._solution_validation_limits(options.tolerances,
                      Dict(:r_D_candidate=>0.0,:r_A_candidate=>0.0)),))
    end # setprecision
end
end
