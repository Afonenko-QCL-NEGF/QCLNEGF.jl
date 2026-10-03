module FrozenScientificFixtures
include("../support/common.jl")

# One immutable small geometry. These operator checks do not run SCBA/Poisson
# and do not establish uniqueness or validate a stationary transport solution.
function fixture()
    p = reference_parameters(Tᴸ = 70u"K", V_period = 0u"mV")
    n = NumericalParameters(
        N_z=18, N_b=1, P_basis=1, N_E=17,
        E_min=-0.15u"eV", E_max=0.15u"eV", M_E=0.02u"eV",
        k_max=0.2u"nm^-1", N_k=2, N_φ=8, qz_max=2u"nm^-1",
        N_qz=5, η_seed=2u"meV",
    )
    s=ScaleSystem()
    g=build_grids(p,n,s)
    profiles=build_profiles(p,n,s,g)
    basis=build_localized_basis(p,n,s,g,profiles)
    return p,n,s,g,profiles,basis
end

function components(support)
    labels=zeros(Int,size(support,1))
    count=0
    for seed in eachindex(labels)
        labels[seed]==0 || continue
        count+=1
        pending=[seed]
        labels[seed]=count
        while !isempty(pending)
            source=pop!(pending)
            for destination in eachindex(labels)
                if labels[destination]==0 && (support[destination,source] || support[source,destination])
                    labels[destination]=count
                    push!(pending,destination)
                end
            end
        end
    end
    return labels
end

@testset "Fixed spectral energy support and sector balance" begin
    # Exact synthetic shifts isolate support from microscopic vertex errors.
    # Positive scalar spectra/kernels make the actual contraction support equal
    # to the remap support; this is a declared one-basis frozen operator fixture.
    energy=collect(-8.0:8.0)
    # Endpoint cells end at the first/last node, as specified by the
    # finite-volume remap; their widths are half the interior width.
    weights=vcat(0.5,ones(15),0.5)
    phonon=QCLNEGF.build_conservative_shift_pair(energy,weights,2.0)
    spectral=reshape(ComplexF64.(1 .+ 0.03 .* (1:17)),17,1,1,1)
    kernel=ones(ComplexF64,1,1,1,1,1,1)
    for (bias,expected_sectors) in ((0.0,2),(0.5,1))
        field=QCLNEGF.build_conservative_shift_pair(energy,weights,bias)
        transfer=zeros(17,17)
        for source in 1:17
            occupation=zeros(ComplexF64,17,1,1,1)
            occupation[source]=spectral[source]
            shifted=QCLNEGF.apply_energy_shift(phonon.plus+phonon.minus+field.plus+field.minus,occupation)
            contraction=QCLNEGF._static_contraction(kernel,shifted,[1.0])
            transfer[:,source]=weights .* real.(vec(spectral).*vec(contraction))
        end
        @test minimum(transfer)>=0
        labels=components(transfer .> 0)
        @test length(unique(labels))==expected_sectors
        occupation=collect(range(0.15,0.85;length=17))
        flux=transfer .* reshape(occupation,1,17) .* reshape(1 .-occupation,17,1)
        balance=vec(sum(flux;dims=2))-vec(sum(flux;dims=1))
        # The finite-window exterior is excluded explicitly. This balance is
        # for the represented pair transfer, not a claim about missing tails.
        for sector in unique(labels)
            @test abs(sum(balance[labels .==sector]))<=64eps(Float64)*sum(abs,flux)
        end
        println("[frozen-support] bias_shift=",bias," sectors=",expected_sectors,
                " edges=",count(>(0),transfer)," total_balance=",sum(balance))
    end
end

@testset "Frozen impurity occupation-weighted contraction" begin
    p,n,s,g,profiles,basis=fixture()
    adaptive,diagnostic=QCLNEGF.build_kernels_exact_parallel(
        p,n,ScatteringOptions(LO=false,acoustic=false,impurity=true,IFR=false,alloy=false),
        s,g,profiles,basis; impurity_angular=:adaptive,angular_tolerance=1e-5,
    )
    raw=adaptive.qᴷ[:impurity].*adaptive.K[:impurity]
    occupation=zeros(ComplexF64,17,2,1,1)
    for e in 1:17,m in 1:2
        occupation[e,m,1,1]=im*(0.1+0.02e)*(1+0.1m)
    end
    target=QCLNEGF._static_contraction(raw,occupation,g.wᵏ)
    oracle=Array{ComplexF64,4}[]
    for nodes in (2048,4096)
        angles=[2π*(j-0.5)/nodes for j in 1:nodes]
        positive_grids=ModelGrids((name===:φ ? angles : getfield(g,name) for name in fieldnames(ModelGrids))...)
        # Literal donor/vertex loops + positive midpoint quadrature are
        # independent of the optimized evaluator and adaptive Simpson rule.
        literal=QCLNEGF.impurity_kernel(p,s,positive_grids,profiles,basis;angular=:legacy_uniform)
        push!(oracle,QCLNEGF._static_contraction(literal,occupation,g.wᵏ))
    end
    scale=norm(last(oracle))
    @test scale>0
    quadrature_change=norm(oracle[2]-oracle[1])/scale
    contraction_error=norm(target-oracle[2])/scale
    @test quadrature_change<1e-5
    @test contraction_error<5e-5
    @test diagnostic.mechanisms[:impurity].angular_quadrature_checked
    println("[frozen-impurity] oracle_change=",quadrature_change,
            " weighted_contraction_error=",contraction_error," scale=",scale)
end

@testset "Finite window PV has an independent analytic element oracle" begin
    gamma=reshape(ComplexF64[0,1,0],3,1,1,1)
    result=QCLNEGF.product_integration_hilbert_transform(gamma,[0.0,1.0,2.0],[0.5,1.0,0.5])
    @test real.(vec(result)) ≈ [-log(2)/π,0,log(2)/π] atol=2e-15
    println("[frozen-PV] analytic_triangle_error=",maximum(abs.(real.(vec(result)).-[-log(2)/π,0,log(2)/π])))
end
end
