"""
    retarget_problem(problem; V_period=nothing, F_bias=nothing,
                     Tᴸ=problem.physical.Tᴸ,
                     Tᴸᴼ=problem.physical.Tᴸᴼ,
                     validate_static=true, energy_shift=nothing)

Create an operating-point copy of an already built problem while reusing the
spatial grid, material profiles, localized basis and normalized microscopic
kernels.  Exactly one of `V_period` and `F_bias` may be supplied; omitting both
keeps the field. The two period-shift operators are rebuilt because the
field-periodic energy shift changes. The original energy discretization is
preserved unless `energy_shift` explicitly selects another one. The LO matrices
are reused when that discretization is unchanged; otherwise both field and LO
operators are rebuilt together.

For the implemented elastic-equipartition acoustic channel the normalized
kernel is temperature independent and its exact scale is proportional to
`Tᴸ`; consequently only `qᴷ[:acoustic]` is rescaled.  LO temperature enters
the Bose factor during every SCBA candidate and requires no kernel rebuild.
No other material or scattering parameter is silently changed.

See [Field-periodic closure](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/main/docs/src/theory/06_periodicity.md),
[Microscopic kernels](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/main/docs/src/theory/08_kernels.md), and
[production sweep construction](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/main/docs/src/theory/19_production.md).
"""
function retarget_problem(
    problem::NEGFProblem;
    V_period = nothing,
    F_bias = nothing,
    Tᴸ = problem.physical.Tᴸ,
    Tᴸᴼ = problem.physical.Tᴸᴼ,
    validate_static::Bool = true,
    energy_shift::Union{Nothing,Symbol} = nothing,
)
    energy_shift === nothing ||
        energy_shift in (:dense, :sparse_plan, :conservative_pair) ||
        throw(
            ArgumentError(
                "energy_shift must be :dense, :sparse_plan, or :conservative_pair",
            ),
        )
    discretization =
        energy_shift === nothing ? problem.energy_shift_discretization :
        energy_shift === :conservative_pair ? :finite_volume_piecewise_constant :
        :nodal_linear
    V_period !== nothing &&
        F_bias !== nothing &&
        throw(ArgumentError("supply either V_period or F_bias, not both"))
    Tnew = _temperature(Tᴸ)
    TLOnew = _temperature(Tᴸᴼ)
    Tnew > 0u"K" && TLOnew > 0u"K" ||
        throw(ArgumentError("lattice and LO temperatures must be positive"))

    p = problem.physical
    Lp = period_length(p)
    field = if F_bias !== nothing
        _field(F_bias)
    elseif V_period !== nothing
        _field(uconvert(u"V/m", V_period / Lp))
    else
        p.F_bias
    end
    field ≥ 0u"V/m" || throw(ArgumentError("reference design sign convention requires F_bias ≥ 0"))

    pnew = PhysicalParameters(
        copy(p.layers),
        p.N_dop²ᴰ,
        p.f_ion,
        field,
        p.z₀,
        p.E_ref,
        Tnew,
        TLOnew,
        p.ε_s,
        p.ε_∞,
        p.ħωᴸᴼ,
        p.q_s,
        p.qᴸᴼ_s,
        p.Δᴵᶠᴿ,
        p.Λᴵᶠᴿ,
        p.Ξ,
        p.ρ_m,
        p.v_s,
        p.ΔV_alloy,
        p.Ω₀,
        p.g_s,
        copy(p.interfaces),
    )

    qᴷ = copy(problem.kernels.qᴷ)
    if :acoustic in problem.kernels.enabled
        qᴷ[:acoustic] *= _kelvin(Tnew) / _kelvin(p.Tᴸ)
    end
    kernels =
        KernelSet(problem.kernels.K, qᴷ, problem.kernels.Fᴸᴼ, copy(problem.kernels.enabled))
    sp = _scaled_physics(pnew, problem.scales)
    ε = problem.grids.ε
    compact =
        energy_shift === :sparse_plan ||
        (energy_shift === nothing && problem.W₊ᴱᵖ isa EnergyShiftPlan)
    shift_builder = compact ? build_shift_plan : build_shift_matrix
    if discretization === :finite_volume_piecewise_constant
        field_shift = build_conservative_shift_pair(ε, problem.grids.wᴱ, sp.Eᵖ)
        Wplus, Wminus = field_shift.plus, field_shift.minus
    else
        Wplus, Wminus = shift_builder(ε, +sp.Eᵖ), shift_builder(ε, -sp.Eᵖ)
    end
    WLOplus, WLOminus =
        if discretization === problem.energy_shift_discretization &&
           compact == (problem.W₊ᴸᴼ isa EnergyShiftPlan)
            problem.W₊ᴸᴼ, problem.W₋ᴸᴼ
        elseif discretization === :finite_volume_piecewise_constant
            phonon_shift = build_conservative_shift_pair(ε, problem.grids.wᴱ, sp.ħωᴸᴼ)
            phonon_shift.plus, phonon_shift.minus
        else
            shift_builder(ε, +sp.ħωᴸᴼ), shift_builder(ε, -sp.ħωᴸᴼ)
        end
    retargeted = NEGFProblem(
        pnew,
        problem.numerical,
        problem.scattering,
        problem.scales,
        problem.grids,
        problem.profiles,
        problem.basis,
        kernels,
        Wplus,
        Wminus,
        WLOplus,
        WLOminus,
        discretization,
        problem.models,
    )
    if validate_static
        report = validate_problem(retargeted)
        report.passed || throw(
            ArgumentError(
                "retargeted model validation failed: " * join(report.messages, "; "),
            ),
        )
    end
    return retargeted
end

"""
    retarget_production_cache(cache, problem)

Reuse BLAS-ready kernel operators at a new bias/temperature operating point
and rebuild only the four `O(N_E)` shift plans.  This is exact provided
[`retarget_problem`](@ref) created `problem`, so the normalized six-index
kernels and all array dimensions are unchanged.

See [Production sweep, warm start, and recovery](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/main/docs/src/theory/19_production.md).
"""
function retarget_production_cache(cache::ProductionCache, problem::NEGFProblem)
    _cache_source_contract(cache, problem)
    Nk, Nb = problem.numerical.N_k, problem.numerical.N_b
    for (name, operator) in cache.kernels
        operator.N_k == Nk && operator.N_b == Nb ||
            throw(DimensionMismatch("cached kernel $name has wrong dimensions"))
    end
    sp = _scaled_physics(problem.physical, problem.scales)
    ε = problem.grids.ε
    return ProductionCache(
        build_shift_plan(ε, +sp.Eᵖ),
        build_shift_plan(ε, -sp.Eᵖ),
        build_shift_plan(ε, +sp.ħωᴸᴼ),
        build_shift_plan(ε, -sp.ħωᴸᴼ),
        cache.kernels,
        cache.fft_hilbert_plan,
        cache.product_hilbert_operator,
        cache.estimate,
        cache.source_grids,
        cache.source_kernel_arrays,
        cache.algorithms,
    )
end

