module IndependentNonlinearAndContinuumAnchors
include("../support/common.jl")
const BN = QCLNEGF

numerical_variant(n; changes...) = NumericalParameters(;
    (
        name => get(changes, name, getfield(n, name)) for
        name in fieldnames(NumericalParameters)
    )...,
)

function scalar_family(value)
    # Equilibrium Keldysh partners keep the production mixer in its physical cone.
    gamma = -2imag(value)
    return SelfEnergyFamily(
        fill(ComplexF64(value), 1, 1, 1, 1),
        fill(0.3im*gamma, 1, 1, 1, 1),
        fill(-0.7im*gamma, 1, 1, 1, 1),
    )
end

function nonlinear_scba(z, lambda, seed; mixing = :linear, maxiter = 1000)
    state = [scalar_family(-im*seed)]
    workspace = BN._AndersonWorkspace()
    algorithms = AlgorithmOptions(mixing = mixing, anderson_history_depth = 4)
    production = ProductionOptions(algorithms = algorithms, worker_count = 1)
    tolerance=SolverTolerances(r_Σ = 1e-10)
    policy=ConvergencePolicy(stagnation_window = 0, stagnation_relative_improvement = 0.0)
    streak=0
    history = Float64[]
    green = zero(ComplexF64)
    for iteration = 1:maxiter
        G, _, _ = BN._retarded_green_production(
            [real(z)],
            zeros(ComplexF64, 1, 1, 1),
            state[1].Σᴿ;
            η = imag(z),
            options = production,
        )
        green = only(G)
        candidate = [scalar_family(lambda*green)]
        raw = BN._selfenergy_residual(state[1], candidate[1])
        push!(history, raw)
        # This oracle harness exercises the actual production Dyson, candidate
        # family, raw residual and mixer; analytic truth is the quadratic root.
        row=SCBAIteration(
            iteration,
            0.0,
            0.0,
            0.0,
            raw,
            0.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            NaN,
            NaN,
            NaN,
            NaN,
            nothing,
            nothing,
        )
        assessment=scba_convergence_assessment(row, tolerance, policy)
        streak=assessment.passed ? streak+1 : 0
        if streak>=policy.required_consecutive_scba_passes
            return (; green, history, converged = true)
        end
        if mixing === :anderson
            BN._anderson_mix!(workspace, state, candidate, algorithms, 0.6, production)
            @test length(workspace.states)<=4
        else
            BN._mix_family_production!(state[1], candidate[1], 0.6, production)
        end
        @test BN._anderson_causal(state)
    end
    return (; green, history, converged = false)
end

@testset "Nonlinear SCBA causal root and honest finite budget" begin
    z = 0.2+0.3im
    lambda=0.04
    roots = ((z-sqrt(z*z-4lambda))/(2lambda), (z+sqrt(z*z-4lambda))/(2lambda))
    oracle = only(filter(root->imag(root)<0, roots))
    for method in (:linear, :anderson), seed in (0.01, 0.1, 0.5)
        result = nonlinear_scba(z, lambda, seed; mixing = method)
        @test result.converged
        @test result.green ≈ oracle rtol=2e-9
        @test imag(result.green)<0
        @info "Nonlinear analytic anchor" method seed candidate_calls=length(result.history)
    end
    bounded = nonlinear_scba(z, lambda, 0.5; maxiter = 1)
    @test !bounded.converged
    @test length(bounded.history)==1
end

@testset "Production Lorentzian quadrature separates window and phase errors" begin
    for phase in (0.0, 0.5)
        errors=Float64[]
        for nodes in (101, 201, 401)
            lower, upper=-0.6, 0.8
            energy=collect(range(lower, upper; length = nodes))
            step=(upper-lower)/(nodes-1)
            centre=0.1+phase*step
            width=0.035
            G, _, _=BN._retarded_green_production(
                energy,
                reshape(ComplexF64[centre], 1, 1, 1),
                fill(-im*width, nodes, 1, 1, 1);
                options = ProductionOptions(worker_count = 1),
            )
            weight=fill(step, nodes)
            weight[[1, end]]./=2
            integral=dot(weight, real.(spectral_function(G)[:, 1, 1, 1]))/(2π)
            exact=(atan((upper-centre)/width)-atan((lower-centre)/width))/π
            push!(errors, abs(integral-exact))
            @test 0 < integral < 1
        end
        @test errors[3] < errors[2] < errors[1]
        @test errors[2]/errors[3] > 3.5
    end
end

@testset "Annular momentum weights converge to the finite 2D Fermi integral" begin
    p=reference_parameters(Tᴸ = 70u"K")
    s=ScaleSystem()
    kBT=BN._electronvolts(CODATA.kᴮₑᵥ*p.Tᴸ)/s.E₀_eV
    mass=0.067*CODATA.m₀
    a=BN._electronvolts(CODATA.ħ^2*(inv(s.L₀_m)*u"m^-1")^2/(2mass))/s.E₀_eV
    mu=0.02/s.E₀_eV
    errors=Float64[]
    for nodes in (65, 129, 257)
        n=numerical_variant(tutorial_numerics(); N_k = nodes, k_max = 0.25u"nm^-1")
        grids=build_grids(p, n, s)
        occupations=[BN._fermi(a*k*k, mu, kBT) for k in grids.κ]
        numerical=p.g_s*dot(grids.wᵏ, occupations)
        exact=p.g_s*kBT/(4π*a)*(log1p(exp(mu/kBT))-log1p(exp((mu-a*last(grids.κ)^2)/kBT)))
        tail=p.g_s*kBT/(4π*a)*log1p(exp((mu-a*last(grids.κ)^2)/kBT))
        @test tail>0
        @test exact+tail ≈ p.g_s*kBT/(4π*a)*log1p(exp(mu/kBT))
        push!(errors, abs(numerical-exact)/exact)
    end
    @test errors[3]<1e-4
    @test errors[1]/errors[2]>3.5
    @test errors[2]/errors[3]>3.5
end

@testset "Continuum sinusoidal Poisson source has second-order spatial convergence" begin
    p=reference_parameters(Tᴸ = 70u"K")
    s=ScaleSystem()
    errors=Float64[]
    for nodes in (32, 64, 128)
        n=numerical_variant(tutorial_numerics(); N_z = nodes)
        grids=build_grids(p, n, s)
        L=sum(grids.wˣ)
        epsilon=12.9
        exact=0.01 .* sin.(2π .* grids.x ./ L)
        donors=fill(1e3, nodes)
        profiles=MaterialProfiles(
            zeros(nodes),
            fill(0.067, nodes),
            fill(0.067, nodes),
            fill(epsilon, nodes),
            donors,
            zeros(nodes),
            ones(Int, nodes),
        )
        # Independently differentiate the continuum function. Never use L*U
        # from the discretized matrix to manufacture this source.
        density=donors .+ epsilon*(2π/L)^2/s.λ_P .* exact
        solution, zeta, residual, neutrality=solve_periodic_poisson(
            profiles,
            grids,
            s,
            density,
        )
        push!(errors, norm(solution-exact)/norm(exact))
        @test abs(sum(solution))<1e-11
        @test abs(zeta)<1e-10
        @test residual<1e-9
        @test neutrality<1e-12
    end
    @test errors[1]/errors[2]>3.9
    @test errors[2]/errors[3]>3.9
end
end
