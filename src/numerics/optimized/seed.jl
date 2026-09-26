@doc raw"""
    _seed_number_weights_production(Gᴿ, grids, g_s; thread_trace=nothing)

Precompute the exact discrete sheet-number contribution of every energy node
for the equilibrium seed.  If ``A=i(Gᴿ-Gᴬ)``, the returned vector is

```math
q_e = \frac{g_s}{2\pi}\,w^E_e
      \sum_m w^k_m\,\operatorname{tr} A_{em},
```

so that the seed number functional is the scalar sum
``N(\mu)=\sum_e q_e f(E_e;\mu)``.  This is the same finite quadrature used by
[`number_functional`](@ref), with the energy-only Fermi factor contracted
after the momentum and basis traces.  It avoids constructing an
`(N_E,N_k,N_b,N_b)` lesser array for every chemical-potential trial.

`Gᴿ` is a dense complex array of shape `(N_E,N_k,N_b,N_b)`; the result is a
real dimensionless vector of length `N_E`.  Independent energy rows are
threaded.  The internal `thread_trace` hook, when supplied, must have length
`N_E` and records the Julia thread used for each row; it exists only for
regression testing of the production schedule.

See [Green functions and density](@ref theory-greens), especially
[the fixed-density normalization](@ref eq-number-normalization), and the
[production seed](@ref theory-production).
"""
function _seed_number_weights_production(
    Gᴿ::AbstractArray{<:Number,4},
    grids::ModelGrids,
    g_s::Integer;
    options::ProductionOptions = ProductionOptions(),
    thread_trace::Union{Nothing,AbstractVector{<:Integer}} = nothing,
)
    NE, Nk, Nb, Nb2 = size(Gᴿ)
    Nb == Nb2 || throw(DimensionMismatch("Gᴿ blocks must be square"))
    length(grids.ε) == NE && length(grids.wᴱ) == NE ||
        throw(DimensionMismatch("energy grid differs from Gᴿ"))
    length(grids.κ) == Nk && length(grids.wᵏ) == Nk ||
        throw(DimensionMismatch("momentum grid differs from Gᴿ"))
    g_s > 0 || throw(ArgumentError("g_s must be positive"))
    thread_trace === nothing ||
        length(thread_trace) == NE ||
        throw(DimensionMismatch("thread_trace must have length N_E"))

    q = zeros(Float64, NE)
    prefactor = Float64(g_s) / (2π)
    _production_parallel_linear!(NE, options) do e
        # A diagonal element is evaluated with exactly the defining
        # i(Gᴿ-Gᴬ) expression.  Its imaginary roundoff is identically zero,
        # while taking `real` makes the real-valued storage contract explicit.
        local momentum_sum = 0.0
        @inbounds for m = 1:Nk
            local trace_A = 0.0
            for a = 1:Nb
                local z = Gᴿ[e, m, a, a]
                trace_A += real(im * (z - conj(z)))
            end
            momentum_sum += grids.wᵏ[m] * trace_A
        end
        q[e] = prefactor * grids.wᴱ[e] * momentum_sum
        thread_trace === nothing || (thread_trace[e] = Base.Threads.threadid())
    end
    all(isfinite, q) || throw(DomainError(q, "seed spectral number weights must be finite"))
    return q
end

function _seed_number_from_weights(
    μ::Real,
    ε::AbstractVector{<:Real},
    q::AbstractVector{<:Real},
    kBT::Real,
)
    length(ε) == length(q) ||
        throw(DimensionMismatch("energy nodes and seed weights differ"))
    kBT > 0 || throw(ArgumentError("kBT must be positive"))
    value = 0.0
    @inbounds for e in eachindex(ε, q)
        value += q[e] * _fermi(ε[e], μ, kBT)
    end
    return value
end

@doc raw"""
    _seed_lesser_production!(Gˡ, Gᴿ, ε, μ, kBT;
                             thread_trace=nothing)

Fill the final equilibrium lesser seed in place from
``G^<(E_e,k_m)=i f(E_e;\mu) A(E_e,k_m)`` and
``A=i(G^R-G^A)``.  Blocks `(e,m)` are independent and are distributed over
Julia threads.  All arrays are dimensionless and have shape
`(N_E,N_k,N_b,N_b)`.

The optional internal `thread_trace` vector has length `N_E*N_k` and records
the thread that wrote each block for a scheduling regression test.

See [Green functions and density](@ref theory-greens) and the
[production backend](@ref theory-production).
"""
function _seed_lesser_production!(
    Gˡ::AbstractArray{ComplexF64,4},
    Gᴿ::AbstractArray{<:Number,4},
    ε::AbstractVector{<:Real},
    μ::Real,
    kBT::Real;
    options::ProductionOptions = ProductionOptions(),
    thread_trace::Union{Nothing,AbstractVector{<:Integer}} = nothing,
)
    size(Gˡ) == size(Gᴿ) || throw(DimensionMismatch("Gˡ and Gᴿ arrays differ"))
    NE, Nk, Nb, Nb2 = size(Gᴿ)
    Nb == Nb2 || throw(DimensionMismatch("Gᴿ blocks must be square"))
    length(ε) == NE || throw(DimensionMismatch("energy axis differs"))
    kBT > 0 || throw(ArgumentError("kBT must be positive"))
    thread_trace === nothing ||
        length(thread_trace) == NE * Nk ||
        throw(DimensionMismatch("thread_trace must have length N_E*N_k"))

    _production_parallel_linear!(NE * Nk, options) do linear
        local e = (linear - 1) % NE + 1
        local m = (linear - 1) ÷ NE + 1
        local f = _fermi(ε[e], μ, kBT)
        @inbounds for b = 1:Nb, a = 1:Nb
            # i*A = -(Gᴿ-Gᴬ), evaluated directly without a temporary A.
            Gˡ[e, m, a, b] = -f * (Gᴿ[e, m, a, b] - conj(Gᴿ[e, m, b, a]))
        end
        thread_trace === nothing || (thread_trace[linear] = Base.Threads.threadid())
    end
    return Gˡ
end

"""
    _seed_green_production(problem, h,
                           options=ProductionOptions(); thread_trace=nothing)

Construct the fixed-density equilibrium Green-function seed for the
production backend.  The threaded production Dyson primitive is evaluated
once.  A real vector of `N_E` exact discrete spectral weights then replaces
the former full `Gˡ` allocation at every chemical-potential bisection step.
Only the final `Gˡ` array is constructed, in parallel over `(E_e,k_m)`.

The result `(green, μ₀)` has the same physical and array contract as the
educational `_seed_green`: `green` is a [`GreenState`](@ref) whose four complex
fields have shape `(N_E,N_k,N_b,N_b)`, and `μ₀` is a dimensionless energy.
No momentum averaging or approximation of the discrete quadrature is made.

For tests, `thread_trace=(weights=..., lesser=...)` may be supplied with
integer vectors of lengths `N_E` and `N_E*N_k`; production callers omit it.

See [the SCBA fixed point](@ref theory-scba),
[Green functions and density](@ref theory-greens), and the
[production backend](@ref theory-production).
"""
function _seed_green_production(
    problem::NEGFProblem,
    h::AbstractArray{<:Number,3},
    options::ProductionOptions = ProductionOptions();
    thread_trace = nothing,
)
    n = problem.numerical
    shape = (n.N_E, n.N_k, n.N_b, n.N_b)
    size(h) == (n.N_k, n.N_b, n.N_b) || throw(DimensionMismatch("h has wrong shape"))
    _check_production_options(options)
    weights_trace = thread_trace === nothing ? nothing : thread_trace.weights
    lesser_trace = thread_trace === nothing ? nothing : thread_trace.lesser

    zeroΣ = zeros(ComplexF64, shape)
    η = effective_seed_parameters(problem).effective_scaled
    Gᴿ, κD, scaleD = _retarded_green_production(problem.grids.ε, h, zeroΣ; η = η, options)
    q = _seed_number_weights_production(
        Gᴿ,
        problem.grids,
        problem.physical.g_s;
        options,
        thread_trace = weights_trace,
    )

    kBT = _electronvolts(CODATA.kᴮₑᵥ * problem.physical.Tᴸ) / problem.scales.E₀_eV
    target = _scaled_physics(problem.physical, problem.scales).Nᴰ²ᴰ
    number_at(μ) = _seed_number_from_weights(μ, problem.grids.ε, q, kBT) - target
    μlo = problem.grids.ε[1] - 50kBT
    μhi = problem.grids.ε[end] + 50kBT
    flo, fhi = number_at(μlo), number_at(μhi)
    flo * fhi ≤ 0 || throw(ArgumentError("seed chemical-potential root is not bracketed"))
    μ₀ = find_zero(number_at, (μlo, μhi), Bisection())

    Gˡ = zeros(ComplexF64, shape)
    _seed_lesser_production!(
        Gˡ,
        Gᴿ,
        problem.grids.ε,
        μ₀,
        kBT;
        options,
        thread_trace = lesser_trace,
    )
    # Bisection enforces the number at machine accuracy. A unilateral rescale
    # after this root would spoil exact Fermi bounds in nearly full tail blocks.
    return _green_state_production(Gᴿ, Gˡ, κD, scaleD, options), μ₀
end
