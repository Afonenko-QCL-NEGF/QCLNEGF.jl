"""Select wider physical energy bounds without changing the original ΔE."""
function _adapted_energy_parameters(
    problem::NEGFProblem,
    coverage::RepresentationCoverage,
    policy::DomainAdaptationPolicy;
    measured_failure::Bool = false,
)
    n = expand_energy_window(problem, coverage; maximum_nodes = policy.maximum_energy_nodes)
    measured_failure || return n
    width = problem.numerical.E_max-problem.numerical.E_min
    lower = min(n.E_min, problem.numerical.E_min-policy.growth_fraction*width/2)
    upper = max(n.E_max, problem.numerical.E_max+policy.growth_fraction*width/2)
    step = (problem.numerical.E_max-problem.numerical.E_min)/(problem.numerical.N_E-1)
    nodes = ceil(Int, ustrip(Unitful.NoUnits, (upper-lower)/step))+1
    nodes <= policy.maximum_energy_nodes || throw(
        ArgumentError(
            "measured spectral tails require more energy nodes than the declared adaptation budget",
        ),
    )
    return NumericalParameters(;
        (
            name=>(
                name===:E_min ? lower :
                name===:E_max ? upper : name===:N_E ? nodes : getfield(n, name)
            ) for name in fieldnames(NumericalParameters)
        )...,
    )
end

"""
    solve_adaptive_production(problem; domain_adaptation, kernel_options, ...)

Execute a declared bounded energy-window expansion at preserved energy step.
Every changed domain builds fresh basis projections, kernels and a fresh cache,
then starts with a cold seed. It is never described as checkpoint continuation.
Stored states are sent to `state_observer` before deciding whether to expand.
The latest result contains every domain's metrics, identities and actual cost.
Unmet momentum/basis resolution is reported; this policy adapts energy only.
"""
function solve_adaptive_production(
    problem::NEGFProblem;
    domain_adaptation::DomainAdaptationPolicy = DomainAdaptationPolicy(),
    options::SolverOptions = reference_production_solver_options(),
    production_options::ProductionOptions = ProductionOptions(),
    kernel_options::ProductionKernelOptions = ProductionKernelOptions(),
    state_observer::Union{Nothing,Function} = nothing,
    checkpoint_sink::Union{Nothing,Function} = nothing,
    history_observer::Union{Nothing,Function} = nothing,
    cache::Union{Nothing,ProductionCache} = nothing,
    initial_Uᴴ::Union{Nothing,AbstractVector{<:Real}} = nothing,
    initial_scba::Union{Nothing,SCBAResult} = nothing,
    consume_initial::Bool = false,
    initial_outer_history::Vector{OuterIteration} = OuterIteration[],
    initial_warnings::Vector{Dict{String,Any}} = Dict{String,Any}[],
    resume_scba::Bool = false,
    initial_adaptation = nothing,
)
    domain_adaptation.mode === :none && return _solve_production(
        problem;
        options,
        production_options,
        state_observer,
        checkpoint_sink,
        cache,
        initial_Uᴴ,
        initial_scba = let state = initial_scba
            initial_scba = nothing
            state
        end,
        consume_initial,
        initial_outer_history,
        initial_warnings,
        resume_scba,
        history_observer,
    )
    revisions =
        initial_adaptation===nothing ? Dict{String,Any}[] :
        Dict{String,Any}[Dict{String,Any}(r) for r in initial_adaptation["revisions"]]
    current = problem
    expansions = initial_adaptation===nothing ? 0 : Int(initial_adaptation["expansions"])
    first_resumed_domain = expansions
    function annotated_sink(sink)
        sink===nothing && return nothing
        return function (state)
            state.observables[:adaptation_checkpoint]=Dict(
                "expansions"=>expansions,
                "revisions"=>deepcopy(revisions),
                "policy"=>_restart_contract_value(domain_adaptation),
            )
            return sink(state)
        end
    end
    latest = nothing
    while true
        estimate = representation_coverage(
            current;
            hartree = expansions==0 && initial_Uᴴ !== nothing ? initial_Uᴴ :
                      zeros(current.numerical.N_z),
        )
        # Estimate-driven preflight expansion consumes the same explicit budget
        # as subsequent measured-tail expansion, without an unnecessary solve.
        if !estimate.covered &&
           expansions < domain_adaptation.maximum_expansions &&
           !(initial_scba!==nothing && expansions==first_resumed_domain)
            numerical = _adapted_energy_parameters(current, estimate, domain_adaptation)
            push!(
                revisions,
                Dict(
                    "revision"=>expansions,
                    "source"=>"preflight_estimate",
                    "discretization_id"=>estimate.discretization_id,
                    "basis_id"=>estimate.basis_id,
                    "E_min_eV"=>_electronvolts(current.numerical.E_min),
                    "E_max_eV"=>_electronvolts(current.numerical.E_max),
                    "NE"=>current.numerical.N_E,
                    "required_min_scaled"=>estimate.required_min,
                    "required_max_scaled"=>estimate.required_max,
                    "reason"=>"Hamiltonian spectrum and shift margin exceed current domain",
                    "seed_transfer"=>"none; no state solved",
                ),
            )
        else
            started = time_ns()
            latest = _solve_production(
                current;
                options,
                production_options,
                history_observer,
                state_observer = annotated_sink(state_observer),
                checkpoint_sink = annotated_sink(checkpoint_sink),
                cache = expansions==first_resumed_domain ? cache : nothing,
                initial_Uᴴ = expansions==first_resumed_domain ? initial_Uᴴ : nothing,
                initial_scba = let state =
                        expansions==first_resumed_domain ? initial_scba : nothing
                    initial_scba = nothing
                    state
                end,
                consume_initial = expansions==first_resumed_domain && consume_initial,
                initial_outer_history = expansions==first_resumed_domain ?
                                        initial_outer_history : OuterIteration[],
                initial_warnings = expansions==first_resumed_domain ? initial_warnings :
                                   Dict{String,Any}[],
                resume_scba = expansions==first_resumed_domain && resume_scba,
            )
            initial_scba = nothing
            measured = measured_representation_coverage(current, latest.scba, latest.Uᴴ)
            estimate = measured.estimate
            energy_bad =
                maximum((
                    measured.tails.tail_low,
                    measured.tails.tail_high,
                    measured.tails.edge_gamma,
                    measured.tails.edge_spectral,
                )) > domain_adaptation.tail_threshold
            momentum_bad = measured.tails.tail_k > domain_adaptation.tail_threshold
            push!(
                revisions,
                Dict(
                    "revision"=>expansions,
                    "source"=>"measured_state",
                    "discretization_id"=>estimate.discretization_id,
                    "basis_id"=>estimate.basis_id,
                    "hartree_id"=>estimate.hartree_id,
                    "E_min_eV"=>_electronvolts(current.numerical.E_min),
                    "E_max_eV"=>_electronvolts(current.numerical.E_max),
                    "NE"=>current.numerical.N_E,
                    "seconds"=>(time_ns()-started)/1e9,
                    "status"=>String(latest.status),
                    "strict_accepted"=>latest.converged,
                    "metrics"=>copy(latest.report.metrics),
                    "current_A_per_m2"=>ustrip(
                        u"A/m^2",
                        latest.observables[:electron_flow_current],
                    ),
                    "momentum_tail_unresolved"=>momentum_bad,
                    "seed_transfer"=>"cold rebuild",
                    "reason"=>energy_bad ?
                              "measured energy tails exceed declared threshold" :
                              "energy coverage measured",
                ),
            )
            if (!energy_bad && estimate.covered) ||
               expansions >= domain_adaptation.maximum_expansions ||
               latest.scba.quality===:invalid
                latest.observables[:domain_adaptation] = Dict(
                    "policy"=>_restart_contract_value(domain_adaptation),
                    "revisions"=>revisions,
                    "expansions"=>expansions,
                    "maximum_solver_executions"=>domain_adaptation.maximum_expansions+1,
                    "energy_coverage_met"=>(!energy_bad && estimate.covered),
                    "transfer"=>"cold_seed; no interpolation or exact restart",
                )
                if energy_bad || !estimate.covered
                    push!(
                        latest.observables[:warnings],
                        Dict{String,Any}(
                            "code"=>"DOMAIN_ADAPTATION_BOUNDED",
                            "scope"=>"point",
                            "message"=>"Declared domain expansion budget ended with unresolved energy coverage.",
                        ),
                    )
                end
                return latest
            end
            numerical = _adapted_energy_parameters(
                current,
                estimate,
                domain_adaptation;
                measured_failure = energy_bad,
            )
        end
        built = build_configured_scattering_problem(
            physical = current.physical,
            numerical = numerical,
            scales = current.scales,
            scattering = current.scattering,
            algorithms = production_options.algorithms,
            kernel_options = kernel_options,
            physical_models = current.models,
        )
        current = built.problem
        expansions += 1
    end
end
