module NonzeroPublicCurrentOracle

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

function analytic_basis_phase(problem)
    Nb=problem.numerical.N_b
    full=ComplexF64[1/2 1/sqrt(2);1/sqrt(2) 0;1/2 -1/sqrt(2)]
    literal=full[:,1:Nb]
    overlaps=[dot(literal[:,a],problem.basis.Φ[:,a]) for a in 1:Nb]
    @assert all(z->isapprox(abs(z),1;atol=ATOL,rtol=RTOL),overlaps)
    phase=Diagonal(overlaps./abs.(overlaps))
    b=problem.basis
    @assert isapprox(b.Φ,literal*phase;atol=ATOL,rtol=RTOL)
    @test b.χ ≈ literal*phase atol=ATOL rtol=RTOL
    @test b.localization == :none
    @test b.centres ≈ fill(1.5,Nb) atol=ATOL rtol=RTOL
    @test b.spreads ≈ [1/sqrt(2),1][1:Nb] atol=ATOL rtol=RTOL
    @test b.overlap_eigenvalues ≈ ones(Nb) atol=ATOL rtol=RTOL
    @test b.κ_overlap ≈ 1 atol=ATOL rtol=RTOL
    energies=[(2-sqrt(2))*hopping,2hopping][1:Nb]
    @test b.H₀ ≈ Matrix(Diagonal(ComplexF64.(energies))) atol=ATOL rtol=RTOL
    @test b.Z ≈ phase'*[1.5 -1/sqrt(2);-1/sqrt(2) 1.5][1:Nb,1:Nb]*phase atol=ATOL rtol=RTOL
    @test b.M⁻¹ ≈ Matrix{ComplexF64}(I,Nb,Nb)./mr atol=ATOL rtol=RTOL
    seam=-hopping.*[1/4 1/(2sqrt(2));-1/(2sqrt(2)) -1/2][1:Nb,1:Nb]
    @test b.T₊ ≈ phase'*seam*phase atol=ATOL rtol=RTOL
    @test b.T₋ ≈ b.T₊' atol=ATOL rtol=RTOL
    @test b.window_eigenvalues ≈ energies atol=ATOL rtol=RTOL
    @test b.Egrid₀ ≈ 2hopping atol=ATOL rtol=RTOL
    @test isfinite(b.r_orth) && isfinite(b.r_eigen)
    @test b.r_translation == b.r_translation_nearest == b.r_translation_two_pairs == 0
    return phase
end

function current_fields(problem, phase)
    shape=(3,2,2,2)
    GR,Gl,Gg,A,SR,Sl,Sg=(zeros(ComplexF64,shape) for _ in 1:7)
    condition=zeros(3,2)
    scale=zeros(3,2)
    Q=ComplexF64[2 im;-im 1]
    B=ComplexF64[1 im;-im 2]
    I2=Matrix{ComplexF64}(I,2,2)
    amplitudes=[1 2;3 4;5 6]
    energies=[(2-sqrt(2))*hopping,2hopping]
    for e in 1:3,m in 1:2
        a=amplitudes[e,m]
        literal=(-im.*(a.*Q+I2)./2,im.*a.*Q,-im.*I2,a.*Q+I2,
            -im.*(B+I2./2)./2,im.*I2./2,-im.*B)
        for (target,block) in zip((GR,Gl,Gg,A,SR,Sl,Sg),literal)
            target[e,m,:,:].=phase'*block*phase
            @test phase*target[e,m,:,:]*phase' ≈ block atol=ATOL rtol=RTOL
        end
        xi=problem.grids.ε[e].-energies.-hopping*problem.grids.κ[m]^2
        x,y=xi
        tau=x^2+y^2+21/8
        det2=(x*y-11/16)^2+((5x+3y)/4)^2
        discriminant=tau^2-4det2
        @assert discriminant >= -(ATOL+RTOL*tau^2)
        q=sqrt(max(0,discriminant))
        condition[e,m]=sqrt((tau+q)/(tau-q))
        scale[e,m]=max(sqrt(x^2+9/16)+1/2,sqrt(y^2+25/16)+1/2)
        # Supplemental checks of derived metadata against its declared matrix,
        # not a claim that manufactured GR solves this Dyson equation.
        D=ComplexF64[x+3im/4 -1/2;1/2 y+5im/4]
        @test condition[e,m] ≈ cond(D) atol=ATOL rtol=RTOL
        @test scale[e,m] ≈ opnorm(D,Inf) atol=ATOL rtol=RTOL
        @test GR[e,m,:,:]-GR[e,m,:,:]' ≈ Gg[e,m,:,:]-Gl[e,m,:,:] atol=ATOL rtol=RTOL
        @test A[e,m,:,:] ≈ im.*(GR[e,m,:,:]-GR[e,m,:,:]') atol=ATOL rtol=RTOL
        @test SR[e,m,:,:]-SR[e,m,:,:]' ≈ Sg[e,m,:,:]-Sl[e,m,:,:] atol=ATOL rtol=RTOL
        @test Gl[e,m,:,:]+Gl[e,m,:,:]' ≈ zeros(2,2) atol=ATOL rtol=RTOL
        @test Gg[e,m,:,:]+Gg[e,m,:,:]' ≈ zeros(2,2) atol=ATOL rtol=RTOL
        for H in (A[e,m,:,:],-im.*Gl[e,m,:,:],im.*Gg[e,m,:,:],
                  -im.*Sl[e,m,:,:],im.*Sg[e,m,:,:])
            @test H ≈ H' atol=ATOL rtol=RTOL
            @test eigmin(Hermitian(H)) >= -ATOL
        end
    end
    checked_plain_fields((GR,Gl,Gg,A,SR,Sl,Sg),condition,scale,shape)
    return GreenState(GR,Gl,Gg,A,condition,scale),SelfEnergyFamily(SR,Sl,Sg)
end

@testset "Nonzero current measure and basis coherence: manufactured functional" begin
    problem=manufactured_problem(3,2,:none)
    phase=analytic_basis_phase(problem)
    green,family=current_fields(problem,phase)
    # Literal scalar expansion tr(B*aQ-I/2)=6a-1 gives rows below.
    resolved_expected=[19,43,67]./(8π^2)
    integrated_expected=43/(4π^2)
    resolved_SI=ustrip.(u"A/m^2/eV",energy_resolved_current(problem,family,green))
    resolved=resolved_SI./(J0ref/E0eV)
    @test resolved ≈ resolved_expected atol=ATOL rtol=RTOL
    physical_dE=E0eV.*[0.5,1.0,0.5]
    integrated=sum(physical_dE.*resolved_SI)/J0ref
    @test integrated ≈ integrated_expected atol=ATOL rtol=RTOL
    flux=boundary_flux(problem,family,green)
    @test eSI*ustrip(u"m^-2*s^-1",flux)/J0ref ≈ integrated_expected atol=ATOL rtol=RTOL
    for wrong in (-resolved,2π.*resolved,
                  resolved.-[5,17,29]./(16π^2),
                  resolved.-3 .* [11,23,35]./(16π^2))
        @test !isapprox(wrong,resolved_expected;atol=ATOL,rtol=RTOL)
    end
    @test !isapprox(sum(resolved),integrated_expected;atol=ATOL,rtol=RTOL)
    @test sum(resolved) ≈ 129/(8π^2) atol=ATOL rtol=RTOL
    @test !isapprox(integrated/E0eV,integrated_expected;atol=ATOL,rtol=RTOL)
    @test !isapprox(-integrated,integrated_expected;atol=ATOL,rtol=RTOL)
    @test !isapprox(2π*integrated,integrated_expected;atol=ATOL,rtol=RTOL)

    bonds_expected=fill(15hopping/(8π^2),2)
    bonds=ustrip.(u"A/m^2",bond_current(problem,green.Gˡ))./J0ref
    @test length(bonds)==2 # internal bonds only; no periodic seam.
    @test bonds ≈ bonds_expected atol=ATOL rtol=RTOL
    diagonal=zeros(ComplexF64,3,2,2,2)
    for e in 1:3,m in 1:2
        a=[1 2;3 4;5 6][e,m]
        diagonal[e,m,:,:].=phase'*ComplexF64[2im*a 0;0 im*a]*phase
    end
    diagonal_bonds=ustrip.(u"A/m^2",bond_current(problem,diagonal))./J0ref
    @test diagonal_bonds ≈ zeros(2) atol=ATOL rtol=RTOL
    @test !isapprox(diagonal_bonds,bonds_expected;atol=ATOL,rtol=RTOL)
    truncated=manufactured_problem(3,1,:none)
    analytic_basis_phase(truncated)
    Gltruncated=reshape(ComplexF64.(2im.*[1 2;3 4;5 6]),3,2,1,1)
    truncated_bonds=ustrip.(u"A/m^2",bond_current(truncated,Gltruncated))./J0ref
    @test truncated_bonds ≈ zeros(2) atol=ATOL rtol=RTOL
    @test !isapprox(truncated_bonds,bonds_expected;atol=ATOL,rtol=RTOL)
    # Different Nb is a different projected manufactured state; no charge
    # renormalization or universal uniform-current requirement is asserted.
end
end
