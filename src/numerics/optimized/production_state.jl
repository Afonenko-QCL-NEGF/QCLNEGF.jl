function _current_green_production(
    problem::NEGFProblem,
    h,
    scattering,
    embedding,
    options::ProductionOptions,
)
    shape = (
        problem.numerical.N_E,
        problem.numerical.N_k,
        problem.numerical.N_b,
        problem.numerical.N_b,
    )
    total =
        _total_family_production(scattering, embedding, problem.kernels.enabled, options)
    Gᴿ, κD, scaleD =
        _retarded_green_production(problem.grids.ε, h, total.Σᴿ; η = 0.0, options)
    Gˡraw = zeros(ComplexF64, shape)
    Gᵍraw = zeros(ComplexF64, shape)
    NE, Nk = shape[1], shape[2]
    Nb = shape[3]
    blocks = NE * Nk
    workers = _production_worker_count(options, blocks)
    _production_worker_ranges!(blocks, workers) do _, first_linear, last_linear
        local GR = zeros(ComplexF64, Nb, Nb)
        local ΣL = similar(GR)
        local temporary = similar(GR)
        local result = similar(GR)
        @inbounds for linear = first_linear:last_linear
            local e = (linear - 1) % NE + 1
            local m = (linear - 1) ÷ NE + 1
            for b = 1:Nb, a = 1:Nb
                GR[a, b] = Gᴿ[e, m, a, b]
                ΣL[a, b] = total.Σˡ[e, m, a, b]
            end
            mul!(temporary, GR, ΣL)
            mul!(result, temporary, adjoint(GR))
            for b = 1:Nb, a = 1:Nb
                Gˡraw[e, m, a, b] = result[a, b]
                ΣL[a, b] = total.Σᵍ[e, m, a, b]
            end
            mul!(temporary, GR, ΣL)
            mul!(result, temporary, adjoint(GR))
            for b = 1:Nb, a = 1:Nb
                Gᵍraw[e, m, a, b] = result[a, b]
            end
        end
    end
    target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    number =
        _number_functional_production(Gˡraw, problem.grids, problem.physical.g_s, options)
    λ = number / target
    isfinite(λ) && λ > 0 || throw(DomainError(λ, "kinetic eigenvalue λ must be positive"))
    if options.algorithms.occupation_normalization === :scalar_lesser
        # Explicit legacy comparator. Its reconstructed G> can violate Pauli
        # positivity and is subject to the unchanged physical acceptance gates.
        _production_parallel_linear!(length(Gˡraw), options) do index
            Gˡraw[index] /= λ
        end
        return _green_state_production(Gᴿ, Gˡraw, κD, scaleD, options), λ, total
    end
    holes =
        -_number_functional_production(Gᵍraw, problem.grids, problem.physical.g_s, options)
    weights = _paired_number_coefficients(number, holes, target)
    a, c = weights.occupied_fraction, weights.empty_fraction
    _production_parallel_linear!(length(Gˡraw), options) do index
        lesser, greater = Gˡraw[index], Gᵍraw[index]
        Gˡraw[index] = a*lesser - c*greater
        Gᵍraw[index] = (1-c)*greater - (1-a)*lesser
    end
    return _green_state_production(Gᴿ, Gˡraw, κD, scaleD, options; greater = Gᵍraw),
    λ,
    total
end
