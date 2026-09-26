"""
    reference_parameters(; f_ion=1, F_bias=nothing, V_period=56u"mV",
                       Tᴸ=200u"K", ΔV_alloy=nothing, Ω₀=nothing)

Published reference design layer sequence together with the explicitly chosen material
and scattering parameters of the educational model.  All dimensional inputs
are Unitful quantities; incompatible units raise `Unitful.DimensionError`.
By default the exact `V_period=56 mV` is primary and `F_bias=V_period/Lp`.
An explicit `F_bias` overrides `V_period`.
The optional alloy parameters `ΔV_alloy` and `Ω₀` must either both be omitted
or both be supplied with energy and volume dimensions, respectively.

See [the reference design structure passport](@ref theory-reference2019),
[model conventions](@ref theory-model), and
[the complete input contract](@ref full-input-contract).
"""
function reference_parameters(;
    f_ion::Real = 1.0,
    F_bias = nothing,
    V_period = 56.0u"mV",
    Tᴸ = 200.0u"K",
    Tᴸᴼ = Tᴸ,
    Δᴵᶠᴿ = 0.10u"nm",
    Λᴵᶠᴿ = 9.0u"nm",
    q_s = 0.20u"nm^-1",
    qᴸᴼ_s = 0.20u"nm^-1",
    ΔV_alloy = nothing,
    Ω₀ = nothing,
)

    0 < f_ion ≤ 1 ||
        throw(ArgumentError("f_ion must lie in (0,1] for the fixed-density solver"))
    Tᴸ > 0u"K" && Tᴸᴼ > 0u"K" ||
        throw(ArgumentError("lattice and LO temperatures must be positive"))
    q_s > 0u"m^-1" && qᴸᴼ_s > 0u"m^-1" ||
        throw(ArgumentError("screening wave numbers must be positive"))
    Δᴵᶠᴿ ≥ 0u"m" && Λᴵᶠᴿ > 0u"m" || throw(
        ArgumentError("IFR height must be nonnegative and correlation length positive"),
    )
    (ΔV_alloy === nothing) == (Ω₀ === nothing) ||
        throw(ArgumentError("ΔV_alloy and Ω₀ must be supplied together"))
    if ΔV_alloy !== nothing
        ΔV_alloy > 0u"eV" || throw(ArgumentError("ΔV_alloy must be positive"))
        Ω₀ > 0u"m^3" || throw(ArgumentError("Ω₀ must be positive"))
    end
    barrier =
        (; x_Al = 0.25, Eᶜ = 0.20775u"eV", mᶻᵣ = 0.08775, m_parallelᵣ = 0.08775, εᵣ = 12.9)
    well = (; x_Al = 0.0, Eᶜ = 0.0u"eV", mᶻᵣ = 0.067, m_parallelᵣ = 0.067, εᵣ = 12.9)
    layers = Layer[
        Layer(d = 3.26u"nm", material = :AlGaAs, doped = false; barrier...),
        Layer(d = 7.99u"nm", material = :GaAs, doped = false; well...),
        Layer(d = 1.90u"nm", material = :AlGaAs, doped = false; barrier...),
        Layer(d = 8.40u"nm", material = :GaAs, doped = false; well...),
        Layer(d = 2.90u"nm", material = :GaAs, doped = true; well...),
        Layer(d = 5.16u"nm", material = :GaAs, doped = false; well...),
    ]
    interfaces = _length.([0.0u"nm", 3.26u"nm", 11.25u"nm", 13.15u"nm"])
    Lp = sum(layer.d for layer in layers)
    field = F_bias === nothing ? uconvert(u"V/m", V_period / Lp) : F_bias
    field ≥ 0u"V/m" || throw(ArgumentError("reference design sign convention requires F_bias ≥ 0"))
    return PhysicalParameters(
        layers,
        _sheetdensity(4.5e14u"m^-2"),
        Float64(f_ion),
        _field(field),
        _length(0.0u"nm"),
        _energy(0.0u"eV"),
        _temperature(Tᴸ),
        _temperature(Tᴸᴼ),
        12.9,
        10.89,
        _energy(36.7u"meV"),
        _wavenumber(q_s),
        _wavenumber(qᴸᴼ_s),
        _length(Δᴵᶠᴿ),
        _length(Λᴵᶠᴿ),
        _energy(7.0u"eV"),
        _massdensity(5317.0u"kg/m^3"),
        _speed(5240.0u"m/s"),
        ΔV_alloy === nothing ? nothing : _energy(ΔV_alloy),
        Ω₀ === nothing ? nothing : _volume(Ω₀),
        2,
        interfaces,
    )
end

"""
Reference discretization from the mathematical specification.

See [Indices, grids, and array shapes](@ref array-contracts) and
[the numerical input contract](@ref full-input-contract).
"""
baseline_numerics() = NumericalParameters(
    N_z = 198,
    N_b = 5,
    P_basis = 3,
    E_min = -0.15u"eV",
    E_max = 0.45u"eV",
    N_E = 2401,
    M_E = 56.0u"meV",
    k_max = 0.60u"nm^-1",
    N_k = 61,
    N_φ = 64,
    qz_max = 10.0u"nm^-1",
    N_qz = 401,
    η_seed = 0.5u"meV",
)

"""
Production reference design discretization.  It keeps the reference spatial, energy,
momentum and basis sizes, but reserves an `80 meV` non-periodic shift/trust
margin so the default I--V sweep through `64 mV` per period remains inside
the stored energy window.  The margin changes validation masks only; it does
not change the energy nodes or claim grid convergence.

See [Production backend](@ref theory-production),
[indices, grids, and array shapes](@ref array-contracts), and
[numerical scales and cost](@ref theory-cost).
"""
reference_production_numerics() = NumericalParameters(
    N_z = 198,
    N_b = 5,
    P_basis = 3,
    E_min = -0.15u"eV",
    E_max = 0.45u"eV",
    N_E = 2401,
    M_E = 80.0u"meV",
    k_max = 0.60u"nm^-1",
    N_k = 61,
    N_φ = 64,
    qz_max = 10.0u"nm^-1",
    N_qz = 401,
    η_seed = 0.5u"meV",
)

"""
Small deterministic grid for examples and automated tests.  It changes the discretization,
not the equations, and must never be used for quantitative reference design conclusions.

See [Indices, grids, and array shapes](@ref array-contracts) and the required
[verification and grid-convergence checks](@ref theory-validation).
"""
tutorial_numerics() = NumericalParameters(
    N_z = 48,
    N_b = 3,
    P_basis = 1,
    E_min = -0.10u"eV",
    E_max = 0.28u"eV",
    N_E = 49,
    M_E = 56.0u"meV",
    k_max = 0.45u"nm^-1",
    N_k = 5,
    N_φ = 8,
    qz_max = 6.0u"nm^-1",
    N_qz = 25,
    η_seed = 2.0u"meV",
)

"""
Final acceptance thresholds, direct updates, and normative safety limits.

See [Verification](@ref theory-validation) and
[the fixed-point input contract](@ref full-input-contract).
"""
baseline_options() = SolverOptions()

"""
Damped, strict-tolerance starting point for the production reference design fixed points.
The damping factors are numerical controls, not material parameters, and must
be varied in a convergence study.  No tolerance is relaxed relative to
[`baseline_options`](@ref).

See [the SCBA fixed point](@ref theory-scba),
[the outer Poisson loop](@ref theory-outer-loop), and
[Production backend](@ref theory-production).
"""
reference_production_solver_options() = SolverOptions(
    α_Σ = 0.10,
    α_P = 0.15,
    max_scba = 2000,
    max_poisson = 100,
    tolerances = SolverTolerances(),
    convergence = ConvergencePolicy(),
)

"""
Relaxed mixing and thresholds for algorithmic smoke runs only.

See [the SCBA fixed point](@ref theory-scba) and
[Verification](@ref theory-validation); these relaxed gates are not the
quantitative acceptance contract.
"""
tutorial_options() = SolverOptions(
    α_Σ = 0.35,
    α_P = 0.25,
    max_scba = 80,
    max_poisson = 30,
    convergence = ConvergencePolicy(
        required_consecutive_scba_passes = 1,
        required_consecutive_poisson_passes = 1,
        stagnation_window = 0,
        stagnation_relative_improvement = 0.0,
    ),
    tolerances = SolverTolerances(
        r_D = 1e-10,
        r_A = 1e-6,
        r_K = 5e-4,
        r_Σ = 5e-4,
        r_λ = 5e-4,
        r_P = 1e-10,
        r_U = 5e-4,
        r_n = 5e-4,
        r_neutral = 1e-8,
        r_J = 1e-3,
        r_C = 1e-3,
        r_power = 5e-2,
        r_PSD = 1e-7,
        r_caus = 1e-7,
        r_sum = 5e-2,
        r_obs = 5e-3,
        r_ζ = 1e-7,
        r_imag = 1e-6,
        r_tail = 1.0,
        r_edge = 1.0,
        r_roundoff = 1e-6,
    ),
)

"""
Baseline enabled channels: LO, acoustic, impurity and IFR; alloy off.

See [Microscopic kernels](@ref theory-kernels).
"""
default_scattering() = ScatteringOptions()

period_length(p::PhysicalParameters) = sum(layer.d for layer in p.layers)
ionized_sheet_density(p::PhysicalParameters) = p.f_ion * p.N_dop²ᴰ
