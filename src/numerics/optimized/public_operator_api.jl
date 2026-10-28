"""
    production_retarded_green(energy, h, sigma_R; η=0.0, options)

Evaluate the selected optimized custom inverse, returning all three owned
resolvent/condition/scale arrays. Axes are `(E,k,a,b)` and `(E,k)`; Green units
are inverse input energy units. This operator call is not a scientific solve.
"""
function production_retarded_green(energy::AbstractVector{<:Real},
    h::AbstractArray{<:Number,3}, sigma_R::AbstractArray{<:Number,4};
    η::Real = 0.0, options::ProductionOptions)
    _check_production_options(options)
    return ProductionRetardedGreenEvaluation(
        _retarded_green_production(energy, h, sigma_R; η, options)...)
end

"""
    production_mix_selfenergy!(target, candidate, alpha; options)

Apply the optimized linear family update in place and return the same target.
All three component shapes must agree. The three writable target components
must be pairwise disjoint, and every target/candidate pair must be disjoint
under `Base.mightalias`; refusal occurs before mutation. Candidate is read-only
and may have internal component aliases. No implicit tensor copies are made;
neither SCBA/Anderson acceptance nor a restart receipt is constructed.
"""
function production_mix_selfenergy!(target::SelfEnergyFamily, candidate::SelfEnergyFamily,
    alpha::Real; options::ProductionOptions)
    _check_production_options(options)
    shape = size(target.Σᴿ)
    for field in (:Σᴿ, :Σˡ, :Σᵍ)
        size(getfield(target, field)) == size(getfield(candidate, field)) == shape ||
            throw(DimensionMismatch("self-energy family components differ"))
    end
    for (left, right) in ((:Σᴿ, :Σˡ), (:Σᴿ, :Σᵍ), (:Σˡ, :Σᵍ))
        Base.mightalias(getfield(target, left), getfield(target, right)) &&
            throw(ArgumentError("target self-energy components must be pairwise disjoint"))
    end
    for left in (:Σᴿ, :Σˡ, :Σᵍ), right in (:Σᴿ, :Σˡ, :Σᵍ)
        Base.mightalias(getfield(target, left), getfield(candidate, right)) &&
            throw(ArgumentError("target and candidate must be disjoint"))
    end
    return _mix_family_production!(target, candidate, alpha, options)
end

"""
    production_lo_selfenergies(problem, green, cache; options)

Evaluate the selected optimized LO lesser/greater collision terms with the
existing thermal/nonthermal policy. Validate the actual cache source/shift/
algorithm binding, LO selection and Green axes. Return fresh arrays; cache
scratch may change. No cache build, retarded reconstruction or normalization
occurs, and problem/Green arrays remain read-only.
"""
function production_lo_selfenergies(problem::NEGFProblem, green::GreenState,
    cache::ProductionCache; options::ProductionOptions)
    _check_production_options(options)
    :LO in problem.kernels.enabled || throw(ArgumentError("LO scattering is not enabled"))
    _validated_production_cache(cache, problem, options)
    n = problem.numerical
    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    for field in (:Gᴿ, :Gˡ, :Gᵍ, :A)
        size(getfield(green, field)) == shape || throw(DimensionMismatch("LO Green axes differ"))
    end
    return LOCollisionEvaluation(_production_lo_contraction(cache.kernels[:LO], green,
        problem, cache, problem.kernels.qᴷ[:LO], options)...)
end

"""
    sampled_linewidth_diagnostics(A, total, energy, energy_weights, momentum_weights; policy)

Observe sampled modal Γ/ΔE with the existing spectral midpoint quantile policy.
A/Σ use consistent scaled units, Γ=i(Σᵍ-Σˡ); no hbar or 2π enters this ratio.
Keep all status/reason symbols and unavailable NaNs; a sampled linewidth is not
proof of discretization adequacy. Input arrays are read-only.
"""
function sampled_linewidth_diagnostics(A::AbstractArray{<:Number,4}, total::SelfEnergyFamily,
    energy::AbstractVector{<:Real}, energy_weights::AbstractVector{<:Real},
    momentum_weights::AbstractVector{<:Real}; policy::SCBAPhysicsMarkerPolicy)
    NE, Nk, Nb, Nb2 = size(A)
    Nb == Nb2 && Nb > 0 && NE > 0 && Nk > 0 ||
        throw(DimensionMismatch("linewidth axes must be nonempty with square blocks"))
    length(energy) == length(energy_weights) == NE && length(momentum_weights) == Nk ||
        throw(DimensionMismatch("linewidth quadrature axes differ"))
    for field in (:Σᴿ, :Σˡ, :Σᵍ)
        size(getfield(total, field)) == size(A) ||
            throw(DimensionMismatch("linewidth self-energy axes differ"))
    end
    grids = (; ε = energy, wᴱ = energy_weights, wᵏ = momentum_weights)
    observed = _sampled_linewidth_markers(A, total, grids, policy)
    return SampledLinewidthDiagnostics(observed.q10, observed.q50, observed.q90,
        observed.underresolved, observed.blocks, observed.status, _MARKER_LINEWIDTH_SAMPLING_METHOD)
end

"""
    fermi_occupation(energy, chemical_potential, kBT)

Evaluate the existing saturated Fermi rule for finite real inputs with kBT>0.
All arguments have the same energy coordinates. This new boundary guard does
not change internal `_fermi` IEEE zero-temperature behavior or LO Bose policy.
"""
function fermi_occupation(energy::Real, chemical_potential::Real, kBT::Real)
    all(isfinite, (energy, chemical_potential, kBT)) && kBT > 0 ||
        throw(ArgumentError("Fermi inputs must be finite with positive kBT"))
    return Float64(_fermi(energy, chemical_potential, kBT))
end
