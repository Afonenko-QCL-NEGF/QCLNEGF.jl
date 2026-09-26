"""
    Layer(; d, material, x_Al, Eᶜ, mᶻᵣ, m_parallelᵣ, εᵣ, doped=false)

One epitaxial segment.  `d` and `Eᶜ` are Unitful length and energy
quantities; the constructor converts them to metres and electronvolts.
Relative masses and permittivity are dimensionless.

Physical meaning and the reference design layer ordering are described in
[the reference design structure passport](@ref theory-reference2019); accepted input dimensions
are listed in [the input contract](@ref full-input-contract).
"""
struct Layer
    d::LengthQuantity
    material::Symbol
    x_Al::Float64
    Eᶜ::EnergyQuantity
    mᶻᵣ::Float64
    m_parallelᵣ::Float64
    εᵣ::Float64
    doped::Bool
    function Layer(d, material, x_Al, Eᶜ, mᶻᵣ, m_parallelᵣ, εᵣ, doped)
        isfinite(d) && d>0u"m" ||
            throw(ArgumentError("layer thickness must be finite and positive"))
        isfinite(Eᶜ) || throw(ArgumentError("band edge must be finite"))
        isfinite(x_Al) && 0<=x_Al<=1 ||
            throw(ArgumentError("alloy fraction must lie in [0,1]"))
        all(x->isfinite(x)&&x>zero(x), (mᶻᵣ, m_parallelᵣ, εᵣ)) ||
            throw(ArgumentError("masses and permittivity must be finite and positive"))
        new(
            _length(d),
            material,
            Float64(x_Al),
            _energy(Eᶜ),
            Float64(mᶻᵣ),
            Float64(m_parallelᵣ),
            Float64(εᵣ),
            doped,
        )
    end

end

function Layer(;
    d,
    material::Symbol,
    x_Al::Real,
    Eᶜ,
    mᶻᵣ::Real,
    m_parallelᵣ::Real,
    εᵣ::Real,
    doped::Bool = false,
)
    d > 0u"m" || throw(ArgumentError("layer thickness must be positive"))
    0 ≤ x_Al ≤ 1 || throw(ArgumentError("x_Al must lie in [0,1]"))
    mᶻᵣ > 0 && m_parallelᵣ > 0 || throw(ArgumentError("effective masses must be positive"))
    εᵣ > 0 || throw(ArgumentError("relative permittivity must be positive"))
    return Layer(
        _length(d),
        material,
        Float64(x_Al),
        _energy(Eᶜ),
        Float64(mᶻᵣ),
        Float64(m_parallelᵣ),
        Float64(εᵣ),
        doped,
    )
end

"""
All dimensioned and dimensionless inputs defining one layered transport model.

See [Model and conventions](@ref theory-model),
[reference design structure passport](@ref theory-reference2019), and
[the complete input contract](@ref full-input-contract).
"""
struct PhysicalParameters
    layers::Vector{Layer}
    N_dop²ᴰ::SheetDensityQuantity
    f_ion::Float64
    F_bias::FieldQuantity
    z₀::LengthQuantity
    E_ref::EnergyQuantity
    Tᴸ::TemperatureQuantity
    Tᴸᴼ::TemperatureQuantity
    ε_s::Float64
    ε_∞::Float64
    ħωᴸᴼ::EnergyQuantity
    q_s::WaveNumberQuantity
    qᴸᴼ_s::WaveNumberQuantity
    Δᴵᶠᴿ::LengthQuantity
    Λᴵᶠᴿ::LengthQuantity
    Ξ::EnergyQuantity
    ρ_m::MassDensityQuantity
    v_s::SpeedQuantity
    ΔV_alloy::Union{Nothing,EnergyQuantity}
    Ω₀::Union{Nothing,VolumeQuantity}
    g_s::Int
    interfaces::Vector{LengthQuantity}
    function PhysicalParameters(
        layers,
        N_dop²ᴰ,
        f_ion,
        F_bias,
        z₀,
        E_ref,
        Tᴸ,
        Tᴸᴼ,
        ε_s,
        ε_∞,
        ħωᴸᴼ,
        q_s,
        qᴸᴼ_s,
        Δᴵᶠᴿ,
        Λᴵᶠᴿ,
        Ξ,
        ρ_m,
        v_s,
        ΔV_alloy,
        Ω₀,
        g_s,
        interfaces,
    )
        !isempty(layers) || throw(ArgumentError("physical model needs at least one layer"))
        isfinite(N_dop²ᴰ) && N_dop²ᴰ>0u"m^-2" ||
            throw(ArgumentError("doping must be finite and positive"))
        isfinite(f_ion) && 0<f_ion<=1 ||
            throw(ArgumentError("ionization fraction must lie in (0,1]"))
        isfinite(F_bias) && F_bias>=0u"V/m" ||
            throw(ArgumentError("bias field must be finite and nonnegative"))
        all(isfinite, (z₀, E_ref)) ||
            throw(ArgumentError("energy and position reference must be finite"))
        all(
            x->isfinite(x)&&x>zero(x),
            (Tᴸ, Tᴸᴼ, ε_s, ε_∞, ħωᴸᴼ, q_s, qᴸᴼ_s, Λᴵᶠᴿ, Ξ, ρ_m, v_s),
        ) || throw(ArgumentError("physical scales must be finite and positive"))
        isfinite(Δᴵᶠᴿ) && Δᴵᶠᴿ>=0u"m" ||
            throw(ArgumentError("IFR height must be finite and nonnegative"))
        (ΔV_alloy===nothing)==(Ω₀===nothing) ||
            throw(ArgumentError("alloy energy and volume must be supplied together"))
        ΔV_alloy===nothing ||
            (isfinite(ΔV_alloy)&&ΔV_alloy>0u"eV"&&isfinite(Ω₀)&&Ω₀>0u"m^3") ||
            throw(ArgumentError("alloy parameters must be finite and positive"))
        g_s>=1 || throw(ArgumentError("spin degeneracy must be positive"))
        all(isfinite, interfaces) ||
            throw(ArgumentError("interface positions must be finite"))
        new(
            layers,
            N_dop²ᴰ,
            f_ion,
            F_bias,
            z₀,
            E_ref,
            Tᴸ,
            Tᴸᴼ,
            ε_s,
            ε_∞,
            ħωᴸᴼ,
            q_s,
            qᴸᴼ_s,
            Δᴵᶠᴿ,
            Λᴵᶠᴿ,
            Ξ,
            ρ_m,
            v_s,
            ΔV_alloy,
            Ω₀,
            g_s,
            interfaces,
        )
    end

end

"""
Numerical grid sizes, bounds, trust margin, and seed broadening used by the
direct implementation.

See [Indices, grids, and array shapes](@ref array-contracts) and
[the complete input contract](@ref full-input-contract).
"""
struct NumericalParameters
    N_z::Int
    N_b::Int
    P_basis::Int
    E_min::EnergyQuantity
    E_max::EnergyQuantity
    N_E::Int
    M_E::EnergyQuantity
    k_max::WaveNumberQuantity
    N_k::Int
    N_φ::Int
    qz_max::WaveNumberQuantity
    N_qz::Int
    η_seed::EnergyQuantity
    function NumericalParameters(
        N_z,
        N_b,
        P_basis,
        E_min,
        E_max,
        N_E,
        M_E,
        k_max,
        N_k,
        N_φ,
        qz_max,
        N_qz,
        η_seed,
    )
        N_z>=3 && N_b>=1 && P_basis>=1 || throw(
            ArgumentError("spatial and basis dimensions must be positive and resolved"),
        )
        N_E>=3 && N_k>=2 && N_φ>=4 && N_qz>=3 ||
            throw(ArgumentError("quadrature dimensions are underspecified"))
        all(isfinite, (E_min, E_max, M_E, k_max, qz_max, η_seed)) ||
            throw(ArgumentError("numerical domain values must be finite"))
        E_max>E_min && M_E>0u"eV" && 2M_E<E_max-E_min ||
            throw(ArgumentError("energy window and trusted margin are inconsistent"))
        k_max>0u"m^-1" && qz_max>0u"m^-1" && η_seed>0u"eV" ||
            throw(ArgumentError("wavevector cutoffs and seed broadening must be positive"))
        new(
            N_z,
            N_b,
            P_basis,
            E_min,
            E_max,
            N_E,
            M_E,
            k_max,
            N_k,
            N_φ,
            qz_max,
            N_qz,
            η_seed,
        )
    end

end

function NumericalParameters(;
    N_z::Integer,
    N_b::Integer,
    P_basis::Integer,
    E_min,
    E_max,
    N_E::Integer,
    M_E,
    k_max,
    N_k::Integer,
    N_φ::Integer,
    qz_max,
    N_qz::Integer,
    η_seed,
)
    N_z ≥ 3 || throw(ArgumentError("N_z must be at least 3"))
    N_b ≥ 1 || throw(ArgumentError("N_b must be positive"))
    P_basis ≥ 1 || throw(ArgumentError("P_basis must be positive"))
    N_E ≥ 3 && N_k ≥ 2 && N_φ ≥ 4 && N_qz ≥ 3 ||
        throw(ArgumentError("all quadratures are underspecified"))
    E_max > E_min || throw(ArgumentError("E_max must exceed E_min"))
    M_E > 0u"eV" || throw(ArgumentError("M_E must be positive"))
    2M_E < E_max - E_min ||
        throw(ArgumentError("the Hilbert trust window must be nonempty"))
    k_max > 0u"m^-1" && qz_max > 0u"m^-1" ||
        throw(ArgumentError("momentum cutoffs must be positive"))
    η_seed > 0u"eV" || throw(ArgumentError("η_seed must be positive"))
    return NumericalParameters(
        Int(N_z),
        Int(N_b),
        Int(P_basis),
        _energy(E_min),
        _energy(E_max),
        Int(N_E),
        _energy(M_E),
        _wavenumber(k_max),
        Int(N_k),
        Int(N_φ),
        _wavenumber(qz_max),
        Int(N_qz),
        _energy(η_seed),
    )
end

"""
    SolverTolerances(; ...)

Dimensionless acceptance thresholds for the inner SCBA map, outer Poisson
map, conservation identities, edge tails, and imaginary cancellation tests.
Every enabled criterion must pass; a non-finite residual always fails.
This compositional record is unchecked until the validated `SolverOptions`
constructor verifies every threshold; it is not a standalone solver input.

The individual gates are defined in [Verification](@ref theory-validation)
and tabulated in [the input contract](@ref full-input-contract).
"""
Base.@kwdef struct SolverTolerances
    r_D::Float64 = 1e-11
    r_A::Float64 = 1e-9
    r_K::Float64 = 1e-8
    r_Σ::Float64 = 1e-8
    r_λ::Float64 = 1e-8
    r_P::Float64 = 1e-11
    r_U::Float64 = 1e-7
    r_n::Float64 = 1e-7
    r_neutral::Float64 = 1e-10
    r_J::Float64 = 1e-6
    r_C::Float64 = 1e-8
    r_power::Float64 = 1e-4
    r_PSD::Float64 = 1e-10
    r_caus::Float64 = 1e-10
    r_sum::Float64 = 1e-3
    r_obs::Float64 = 1e-4
    r_ζ::Float64 = 1e-10
    r_imag::Float64 = 1e-10
    r_tail::Float64 = 1e-6
    r_edge::Float64 = 1e-6
    r_roundoff::Float64 = 1e-10
end

"""
    DiagnosticQualityPolicy(; enabled=false, ...)

Explicit approximate acceptance for exploratory studies. The target band is
tried before the coarser band; strict attempts after target acceptance are
bounded. Algebraic gates remain strict and physical quality gates have their
own declared tolerances. Approximate acceptance never sets `converged=true`.
This unchecked compositional record is validated by `SolverOptions` together
with its enclosing convergence policy.
"""
Base.@kwdef struct DiagnosticQualityPolicy
    enabled::Bool = false
    keldysh_threshold::Float64 = 1e-3
    self_energy_threshold::Float64 = 1e-3
    normalization_threshold::Float64 = 1e-4
    required_consecutive_passes::Int = 3
    target_fixed_point_threshold::Float64 = 1e-4
    coarse_wait_iterations::Int = 64
    strict_attempt_iterations::Int = 8
    observable_threshold::Float64 = 1e-3
    positivity_threshold::Float64 = 1e-6
    causality_threshold::Float64 = 1e-10
    outer_threshold::Float64 = 1e-4
end

"""Trend parameters on fresh, unmixed r_K, r_Σ and r_λ histories.

The absolute promising threshold expresses numerical usefulness, never scientific
acceptance. A late plateau is allowed after sufficient reduction. Window means
arithmetic means; regression compares the final window to the best full window.
Contract version 2 applies relative regression only outside the useful band;
an entire window inside that band permits research continuation independently.
"""
Base.@kwdef struct ResearchTrendPolicy
    window::Int = 32
    minimum_reduction::Float64 = 0.5
    promising_threshold::Float64 = 1e-4
    max_regression::Float64 = 2.0
    max_oscillation::Float64 = 10.0

    function ResearchTrendPolicy(
        window::Integer,
        minimum_reduction::Real,
        promising_threshold::Real,
        max_regression::Real,
        max_oscillation::Real,
    )
        window >= 2 || throw(ArgumentError("trend window must be at least two"))
        isfinite(minimum_reduction) && 0 < minimum_reduction < 1 ||
            throw(ArgumentError("minimum_reduction must lie in (0,1)"))
        isfinite(promising_threshold) && promising_threshold > 0 ||
            throw(ArgumentError("promising_threshold must be finite and positive"))
        isfinite(max_regression) && max_regression >= 1 ||
            throw(ArgumentError("max_regression must be finite and at least one"))
        isfinite(max_oscillation) && max_oscillation >= 1 ||
            throw(ArgumentError("max_oscillation must be finite and at least one"))
        new(
            Int(window),
            Float64(minimum_reduction),
            Float64(promising_threshold),
            Float64(max_regression),
            Float64(max_oscillation),
        )
    end
end

"""
    ConvergencePolicy(; ...)

Numerical termination policy shared by the educational and production
solvers.  Residual thresholds remain in [`SolverTolerances`](@ref); this
policy controls only when a sequence of threshold evaluations is accepted
and when a demonstrably stalled SCBA sequence is stopped.
This unchecked compositional record is validated when incorporated into
`SolverOptions`, including direct positional construction of that solver input.

`stagnation_window=0` disables early stagnation termination.  Otherwise the
best normalized fixed-point score (Keldysh, raw self-energy, normalization,
current and population changes) in the newer half of the window must improve
on the older half by at least `stagnation_relative_improvement`. Physical
admissibility gates are diagnosed separately; a saturated positivity defect
cannot hide improvement of the fixed-point iteration.
"""
Base.@kwdef struct ConvergencePolicy
    mode::Symbol = :strict_fail_fast
    trend::ResearchTrendPolicy = ResearchTrendPolicy()
    minimum_scba_iterations::Int = 2
    required_consecutive_scba_passes::Int = 3
    minimum_poisson_iterations::Int = 2
    required_consecutive_poisson_passes::Int = 2
    stagnation_window::Int = 200
    stagnation_relative_improvement::Float64 = 0.01
    diagnostic_quality::DiagnosticQualityPolicy = DiagnosticQualityPolicy()
end

"""
    SolverOptions(; α_Σ=1, α_P=1, max_scba=1000, max_poisson=100, ...)

Controls linear mixing, safety iteration limits, validation tolerances, and
the explicit convergence and energy/momentum tail-window policies.  Mixing
is applied identically to the retarded, lesser, and greater components of
every self-energy family and is recorded in the returned histories.

See [the SCBA fixed point](@ref theory-scba),
[the outer Poisson loop](@ref theory-outer-loop), and
[the input contract](@ref full-input-contract).
"""
Base.@kwdef struct SolverOptions
    α_Σ::Float64 = 1.0
    α_P::Float64 = 1.0
    max_scba::Int = 1000
    max_poisson::Int = 100
    tolerances::SolverTolerances = SolverTolerances()
    convergence::ConvergencePolicy = ConvergencePolicy()
    energy_tail_window_fraction::Float64 = 0.05
    momentum_tail_window_fraction::Float64 = 0.10
    function SolverOptions(
        α_Σ,
        α_P,
        max_scba,
        max_poisson,
        tolerances,
        convergence,
        energy_tail_window_fraction,
        momentum_tail_window_fraction,
    )
        return _check_options(
            new(
                α_Σ,
                α_P,
                max_scba,
                max_poisson,
                tolerances,
                convergence,
                energy_tail_window_fraction,
                momentum_tail_window_fraction,
            ),
        )
    end

end

"""
Boolean switches selecting the five explicitly implemented SCBA channels.

The corresponding microscopic models are defined in
[Microscopic kernels](@ref theory-kernels) and assembled by the
[SCBA map](@ref theory-scba).
"""
Base.@kwdef struct ScatteringOptions
    LO::Bool = true
    acoustic::Bool = true
    impurity::Bool = true
    IFR::Bool = true
    alloy::Bool = false
end

"""
Base units and their derived dimensionless coupling constants.

See [Unitful and the dimensionless core](@ref theory-scaling).
"""
struct ScaleSystem
    E₀::EnergyQuantity
    L₀::LengthQuantity
    m_ref::MassQuantity
    E₀_eV::Float64
    E₀_J::Float64
    L₀_m::Float64
    λ_P::Float64
    J₀::typeof(1.0u"A/m^2")
    function ScaleSystem(E₀, L₀, m_ref, E₀_eV, E₀_J, L₀_m, λ_P, J₀)
        all(x->isfinite(x)&&x>zero(x), (E₀, L₀, m_ref, E₀_eV, E₀_J, L₀_m, λ_P, J₀)) ||
            throw(ArgumentError("solver scales must be finite and positive"))
        isapprox(_electronvolts(E₀), E₀_eV; rtol = 64eps(Float64)) ||
            throw(ArgumentError("E0 eV representation differs"))
        isapprox(ustrip(u"J", uconvert(u"J", E₀)), E₀_J; rtol = 64eps(Float64)) ||
            throw(ArgumentError("E0 joule representation differs"))
        isapprox(_metres(L₀), L₀_m; rtol = 64eps(Float64)) ||
            throw(ArgumentError("L0 metre representation differs"))
        new(E₀, L₀, m_ref, E₀_eV, E₀_J, L₀_m, λ_P, J₀)
    end

end

"""
All quadrature nodes, weights, and the trusted-energy mask in dimensionless
variables.

This is an internal unchecked array representation, exposed for compatibility.
Use the grid builder to construct axes and quadrature together; direct assembly
does not certify correspondence to a `NumericalParameters` value.

Their axes, shapes, and quadrature definitions are specified in
[Indices, grids, and array shapes](@ref array-contracts).
"""
struct ModelGrids
    x::Vector{Float64}
    wˣ::Vector{Float64}
    ε::Vector{Float64}
    wᴱ::Vector{Float64}
    trusted_energy::BitVector
    κ::Vector{Float64}
    wᵏ::Vector{Float64}
    qᶻ::Vector{Float64}
    wᑫᶻ::Vector{Float64}
    φ::Vector{Float64}
end

"""
Piecewise material and ionized-donor profiles on the cell-centred grid.

Internal unchecked array representation. The profile builder owns the grid and
material correspondence; a direct positional call is not a validated input API.

See [the reference design profile](@ref theory-reference2019) and
[the array contract](@ref array-contracts).
"""
struct MaterialProfiles
    Eᶜ::Vector{Float64}
    mᶻᵣ::Vector{Float64}
    m_parallelᵣ::Vector{Float64}
    εᵣ::Vector{Float64}
    Nᴰ::Vector{Float64}
    x_Al::Vector{Float64}
    layer_index::Vector{Int}
end

"""
Fixed cell-local PzP/Löwdin basis and projected one-particle operators.

Internal unchecked numerical representation. Basis builders and static
validation establish dimensions, localization, orthogonality and residuals.

See [BenDaniel–Duke and the fixed PzP basis](@ref theory-basis).
"""
struct BasisData
    localization::Symbol
    Φ::Matrix{ComplexF64}
    χ::Matrix{ComplexF64}
    centres::Vector{Float64}
    spreads::Vector{Float64}
    overlap_eigenvalues::Vector{Float64}
    κ_overlap::Float64
    H₀::Matrix{ComplexF64}
    Z::Matrix{ComplexF64}
    M⁻¹::Matrix{ComplexF64}
    T₊::Matrix{ComplexF64}
    T₋::Matrix{ComplexF64}
    window_eigenvalues::Vector{Float64}
    Egrid₀::Float64
    r_orth::Float64
    r_eigen::Float64
    r_translation::Float64
    r_translation_nearest::Float64
    r_translation_two_pairs::Float64
end

"""
Normalized six-axis kernels in order `(m,m′,a,c,d,b)`.

Internal unchecked array representation. Kernel builders establish the common
grid/basis axes and normalization; the positional container does not certify them.

`K[name]` stores the normalized tensor `Khat_s = Kbar_s/qK_s`, and
`qᴷ[name]` stores the positive scalar `qK_s = maximum(abs, Kbar_s)`.  The
contraction restores this
factor exactly, keeping its summands close to one without changing the map.

See [Microscopic kernels and contractions](@ref theory-kernels) and
[the array contract](@ref array-contracts).
"""
struct KernelSet
    K::Dict{Symbol,Array{ComplexF64,6}}
    qᴷ::Dict{Symbol,Float64}
    Fᴸᴼ::Array{ComplexF64,3}
    enabled::Vector{Symbol}
end

"""Explicit physical approximations; numerical hardware options never select these.

Density-temperature screening uses parabolic 2D compressibility for the donor
potential and homogeneous 3D Debye screening for bulk LO. Kane modifies only
the local transverse kinetic energy. The uniform hot-LO kinetic closure needs
both a measured/declared decay time and density of participating phonon modes.
"""
Base.@kwdef struct PhysicalModelOptions
    screening::Symbol = :fixed
    screening_temperature_K::Float64 = 0.0
    screening_mass_ratio::Float64 = 0.0
    dispersion::Symbol = :parabolic
    nonparabolicity_per_eV::Float64 = 0.0
    lo_population::Symbol = :thermal
    lo_fixed_occupation::Float64 = 0.0
    lo_decay_ps::Float64 = 0.0
    lo_mode_density_per_m3::Float64 = 0.0
    electron_electron::ElectronElectronOptions = ElectronElectronOptions()

    function PhysicalModelOptions(
        screening::Symbol,
        temperature::Real,
        mass::Real,
        dispersion::Symbol,
        alpha::Real,
        lo_population::Symbol,
        occupation::Real,
        decay::Real,
        mode_density::Real,
        electron_electron::ElectronElectronOptions = ElectronElectronOptions(),
    )
        screening in (:fixed, :density_temperature) ||
            throw(ArgumentError("unsupported screening closure"))
        dispersion in (:parabolic, :kane_inplane) ||
            throw(ArgumentError("unsupported dispersion closure"))
        lo_population in (:thermal, :fixed_occupation, :rate_balance) ||
            throw(ArgumentError("unsupported LO population closure"))
        all(
            x -> isfinite(x) && x >= 0,
            (temperature, mass, alpha, occupation, decay, mode_density),
        ) ||
            throw(ArgumentError("physical model parameters must be finite and nonnegative"))
        screening === :fixed ||
            (temperature > 0 && mass > 0) ||
            throw(
                ArgumentError(
                    "density_temperature screening requires electron temperature and band mass",
                ),
            )
        dispersion === :parabolic &&
            alpha != 0 &&
            throw(ArgumentError("nonzero nonparabolicity requires kane_inplane"))
        lo_population === :rate_balance &&
            !(decay > 0 && mode_density > 0) &&
            throw(
                ArgumentError(
                    "LO rate_balance requires decay time and participating mode density",
                ),
            )
        new(
            screening,
            Float64(temperature),
            Float64(mass),
            dispersion,
            Float64(alpha),
            lo_population,
            Float64(occupation),
            Float64(decay),
            Float64(mode_density),
            electron_electron,
        )
    end
end

"""
Complete dimensionless problem after Unitful validation and scaling.

This positional container is an internal unchecked assembled representation,
retained in the facade for existing numerical extensions. Public configured
builders assemble matching axes, units and basis. A hand-built container must
pass static validation before it is treated as a physical problem.

Its construction order is given in [the implementation plan](@ref implementation-plan),
with fields and dimensions in [the input contract](@ref full-input-contract).
"""
struct NEGFProblem{S<:AbstractMatrix{Float64}}
    physical::PhysicalParameters
    numerical::NumericalParameters
    scattering::ScatteringOptions
    scales::ScaleSystem
    grids::ModelGrids
    profiles::MaterialProfiles
    basis::BasisData
    kernels::KernelSet
    W₊ᴱᵖ::S
    W₋ᴱᵖ::S
    W₊ᴸᴼ::S
    W₋ᴸᴼ::S
    # Discretization belongs to the assembled operators, independently of the
    # algorithm used to apply them (dense multiplication or a sparse plan).
    energy_shift_discretization::Symbol
    models::PhysicalModelOptions
end

function NEGFProblem(
    physical,
    numerical,
    scattering,
    scales,
    grids,
    profiles,
    basis,
    kernels,
    W₊ᴱᵖ,
    W₋ᴱᵖ,
    W₊ᴸᴼ,
    W₋ᴸᴼ,
    discretization::Symbol,
)
    return NEGFProblem(
        physical,
        numerical,
        scattering,
        scales,
        grids,
        profiles,
        basis,
        kernels,
        W₊ᴱᵖ,
        W₋ᴱᵖ,
        W₊ᴸᴼ,
        W₋ᴸᴼ,
        discretization,
        PhysicalModelOptions(),
    )
end

# The direct/reference construction contract predates conservative remapping.
# Existing callers supplying the four nodal matrices retain that meaning.
function NEGFProblem(
    physical,
    numerical,
    scattering,
    scales,
    grids,
    profiles,
    basis,
    kernels,
    W₊ᴱᵖ,
    W₋ᴱᵖ,
    W₊ᴸᴼ,
    W₋ᴸᴼ,
)
    return NEGFProblem(
        physical,
        numerical,
        scattering,
        scales,
        grids,
        profiles,
        basis,
        kernels,
        W₊ᴱᵖ,
        W₋ᴱᵖ,
        W₊ᴸᴼ,
        W₋ᴸᴼ,
        :nodal_linear,
    )
end

"""
Dimensionless Green arrays with shape `(N_E,N_k,N_b,N_b)`.

Internal unchecked array representation. Solvers and checkpoint loading bind
these arrays to the actual problem; a positional call alone certifies neither
axis correspondence nor a physical Green state.

See [Green functions and density](@ref theory-greens) and
[the array contract](@ref array-contracts).
"""
struct GreenState
    Gᴿ::Array{ComplexF64,4}
    Gˡ::Array{ComplexF64,4}
    Gᵍ::Array{ComplexF64,4}
    A::Array{ComplexF64,4}
    condition_number::Matrix{Float64}
    dyson_scale::Matrix{Float64}
end

"""
One retarded/lesser/greater self-energy family, all dimensionless.

Internal unchecked mutable array representation. The SCBA map owns its shape
and physical interpretation; direct assembly does not certify a valid state.

See [Microscopic kernels](@ref theory-kernels) and
[the SCBA map](@ref theory-scba).
"""
mutable struct SelfEnergyFamily
    Σᴿ::Array{ComplexF64,4}
    Σˡ::Array{ComplexF64,4}
    Σᵍ::Array{ComplexF64,4}
end

"""Measured PSD witness in solver units, with optional compact reconstruction blocks.

Indices refer to the owning energy/momentum grid. Empty matrices mean that the
historical scalar witness was retained without its reconstruction payload.
Persistence selects first, worst and last payloads; no matrix is reconstructed
from a scalar witness or from a later solver state.
"""
struct SCBAPhysicsWitness
    matrix_kind::Symbol
    energy_index::Int
    momentum_index::Int
    minimum_eigenvalue::Float64
    block_norm::Float64
    backward_error::Float64
    absolute_defect::Float64
    relative_defect::Float64
    ratio::Float64
    hermiticity_defect::Float64
    matrices::Dict{String,Matrix{ComplexF64}}
    eigenvector::Vector{ComplexF64}
end

Base.isequal(a::SCBAPhysicsWitness, b::SCBAPhysicsWitness) = all(
    field->isequal(getfield(a, field), getfield(b, field)),
    fieldnames(SCBAPhysicsWitness),
)

"""Bounded diagnostic work; these settings never change the solver acceptance gates.

Spectral samples are deterministic midpoint quantiles of positive spectral
quadrature weight. Sub-roundoff negative spectral weights are ignored only in
this diagnostic distribution; all unmodified matrices still enter PSD gates.
"""
struct SCBAPhysicsMarkerPolicy
    cadence::Int
    max_spectral_blocks::Int
    relative_mode_weight_floor::Float64
    function SCBAPhysicsMarkerPolicy(;
        cadence::Integer = 10,
        max_spectral_blocks::Integer = 512,
        relative_mode_weight_floor::Real = 1e-10,
    )
        cadence >= 1 || throw(ArgumentError("physics marker cadence must be positive"))
        max_spectral_blocks >= 1 ||
            throw(ArgumentError("spectral sample budget must be positive"))
        0 <= relative_mode_weight_floor < 1 ||
            throw(ArgumentError("relative mode weight floor must lie in [0,1)"))
        new(Int(cadence), Int(max_spectral_blocks), Float64(relative_mode_weight_floor))
    end
end

"""Raw unmixed fixed-point component, using the same Frobenius scale as rΣ."""
struct SCBAChannelMarker
    channel::Symbol
    component::Symbol
    residual_absolute::Float64
    residual_scale::Float64
    residual_relative::Float64
end

"""Collision integral on one fixed normalized Green state, before SCBA mixing.

Signed particle and energy integrals and their absolute incoming/outgoing scales
are dimensionless; the retained imaginary parts expose cancellation/roundoff.
A fresh map and the mixed map have distinct state_kind values.
"""
struct SCBACollisionMarker
    channel::Symbol
    state_kind::Symbol
    particle_signed::Float64
    particle_absolute::Float64
    particle_imaginary::Float64
    energy_signed::Float64
    energy_absolute::Float64
    energy_imaginary::Float64
end

for T in (SCBAChannelMarker, SCBACollisionMarker)
    @eval Base.isequal(a::$T, b::$T) =
        all(field -> isequal(getfield(a, field), getfield(b, field)), fieldnames($T))
end

"""Compact diagnostic measurements of one unmixed Green state, never acceptance gates.

Charges are dimensionless sheet numbers. Corrections and FDT use quadrature-
weighted matrix Frobenius norms. FDT compares both occupied and empty matrices
with a single chemical potential fitted to the *spectral* target number.
Gamma is an energy linewidth, i(SigmaGreater-SigmaLesser), not an angular or
ordinary frequency. Sampled linewidth ratios are diagnostics, not peak FWHM
measurements or a certificate of quadrature convergence.
"""
struct SCBAPhysicalMarkers
    measured_iteration::Int
    raw_hole_charge::Float64
    represented_capacity::Float64
    occupied_fraction_a::Float64
    empty_fraction_c::Float64
    relative_correction_Gn::Float64
    relative_correction_Gp::Float64
    equilibrium_applicable::Bool
    equilibrium_status::Symbol
    equilibrium_mu_eV::Float64
    fdt_raw::Float64
    fdt_normalized::Float64
    equilibrium_abs_current_A_m2::Float64
    delta_energy_eV::Float64
    sampled_gamma_over_dE_q10::Float64
    sampled_gamma_over_dE_q50::Float64
    sampled_gamma_over_dE_q90::Float64
    sampled_spectral_weight_underresolved::Float64
    linewidth_sampled_blocks::Int
    linewidth_status::Symbol
    linewidth_sampling_method::Symbol
    marker_cadence::Int
    linewidth_max_blocks::Int
    relative_mode_weight_floor::Float64
    fresh_map_status::Symbol
    lo_shift_over_dE::Float64
    field_shift_over_dE::Float64
    lo_boundary_occupied_fraction::Float64
    lo_boundary_spectral_fraction::Float64
    field_boundary_occupied_fraction::Float64
    field_boundary_spectral_fraction::Float64
    seed_mu_eV::Float64
    seed_number_ratio::Float64
    fdt_raw_seed_mu::Float64
    fdt_normalized_seed_mu::Float64
    channels::Vector{SCBAChannelMarker}
    collisions::Vector{SCBACollisionMarker}
end

Base.isequal(a::SCBAPhysicalMarkers, b::SCBAPhysicalMarkers) = all(
    field -> isequal(getfield(a, field), getfield(b, field)),
    fieldnames(SCBAPhysicalMarkers),
)

"""
One complete set of residuals for a consistent, pre-mixing SCBA state.

See [the inner SCBA iteration](@ref theory-scba) and
[the validation criteria](@ref theory-validation).
"""
struct SCBAIteration
    ν::Int
    r_D::Float64
    r_A::Float64
    r_K::Float64
    r_Σ::Float64
    r_λ::Float64
    λ::Float64
    r_PSD::Float64
    r_caus::Float64
    r_roundoff::Float64
    r_Jchange::Float64
    r_population::Float64
    J::Float64
    raw_charge::Float64
    target_charge::Float64
    normalized_charge::Float64
    lambda_change::Float64
    witness::Union{Nothing,SCBAPhysicsWitness}
    physical_markers::Union{Nothing,SCBAPhysicalMarkers}

end

"""Algorithm-owned history at an accepted [SCBA boundary](@ref theory-scba), before the next mixing step."""
struct SCBAMixerState
    method::Symbol
    states::Vector{Vector{SelfEnergyFamily}}
    residuals::Vector{Vector{SelfEnergyFamily}}
end
SCBAMixerState() =
    SCBAMixerState(:linear, Vector{SelfEnergyFamily}[], Vector{SelfEnergyFamily}[])

"""Fixed-Hartree [SCBA state](@ref theory-scba), mechanism families, convergence and algorithm history."""
struct SCBAResult
    green::GreenState
    scattering::Dict{Symbol,SelfEnergyFamily}
    embedding::SelfEnergyFamily
    embedding_plus::SelfEnergyFamily
    embedding_minus::SelfEnergyFamily
    history::Vector{SCBAIteration}
    converged::Bool
    status::Symbol
    quality::Symbol
    restart_contract::Union{Nothing,Dict{String,Any}}
    mixer_state::SCBAMixerState

    function SCBAResult(
        green::GreenState,
        scattering::Dict{Symbol,SelfEnergyFamily},
        embedding::SelfEnergyFamily,
        embedding_plus::SelfEnergyFamily,
        embedding_minus::SelfEnergyFamily,
        history::Vector{SCBAIteration},
        converged::Bool,
        status::Symbol,
        quality::Symbol,
        restart_contract::Union{Nothing,Dict{String,Any}} = nothing,
        mixer_state::SCBAMixerState = SCBAMixerState(),
    )
        _check_inner_scba_result_classification(converged, status, quality)
        return new(
            green,
            scattering,
            embedding,
            embedding_plus,
            embedding_minus,
            history,
            converged,
            status,
            quality,
            restart_contract,
            mixer_state,
        )
    end
end

"""Fresh raw self-energy map evaluated on one published Green state."""
struct FixedPointAuditCandidate
    scattering::Dict{Symbol,SelfEnergyFamily}
    embedding::SelfEnergyFamily
    plus::SelfEnergyFamily
    minus::SelfEnergyFamily
    roundoff::Float64
end

"""Particle-number constraint diagnostics in the solver's dimensionless area units."""
struct ChargeConstraintDiagnostics
    raw_charge::Float64
    target_charge::Float64
    normalized_charge::Float64
    lambda::Float64
    lambda_change::Float64
    constraint_residual::Float64
    normalization_residual::Float64
end

"""
One accepted evaluation of the outer Poisson fixed-point map.

See [Periodic Poisson and the outer loop](@ref theory-outer-loop).
"""
struct OuterIteration
    μ::Int
    r_P::Float64
    r_U::Float64
    r_n::Float64
    r_neutral::Float64
    r_J::Float64
    r_Jchange::Float64
    r_population::Float64
    ζ::Float64
    J::Float64
    populations::Vector{Float64}
    density::Vector{Float64}
end

"""
Named validation metrics and all human-readable acceptance failures.

See [Verification](@ref theory-validation).
"""
struct ConvergenceReport
    passed::Bool
    metrics::Dict{Symbol,Float64}
    messages::Vector{String}
end

"""
Complete closed-solver result.

`Uᴴ` and `n` are the dimensionless arrays ``U_H/E_0`` and ``L_0^3 n``;
physical observables are stored in `observables` as Unitful quantities.
`status == :converged` is issued only when both fixed points and the final
independent validation report pass.

See [the closed Poisson–SCBA loop](@ref theory-outer-loop),
[observables and balances](@ref theory-observables), and
[final verification](@ref theory-validation).
"""
struct NEGFSolution
    problem::NEGFProblem
    options::SolverOptions
    Uᴴ::Vector{Float64}
    n::Vector{Float64}
    scba::SCBAResult
    outer_history::Vector{OuterIteration}
    observables::Dict{Symbol,Any}
    report::ConvergenceReport
    converged::Bool
    status::Symbol
end

function Base.show(io::IO, problem::NEGFProblem)
    n = problem.numerical
    print(
        io,
        "NEGFProblem( N_z=$(n.N_z), N_b=$(n.N_b), ",
        "N_E=$(n.N_E), N_k=$(n.N_k), mechanisms=$(problem.kernels.enabled))",
    )
end

function Base.show(io::IO, solution::NEGFSolution)
    print(
        io,
        "NEGFSolution(status=$(solution.status), converged=$(solution.converged), ",
        "outer_iterations=$(length(solution.outer_history)))",
    )
end

"""Explicit bounded domain expansion; each execution cold-rebuilds its operators."""
Base.@kwdef struct DomainAdaptationPolicy
    mode::Symbol = :none
    maximum_expansions::Int = 2
    maximum_energy_nodes::Int = 100_001
    growth_fraction::Float64 = 0.5
    tail_threshold::Float64 = 1e-6
    function DomainAdaptationPolicy(
        mode::Symbol,
        expansions::Integer,
        nodes::Integer,
        growth::Real,
        tail::Real,
    )
        mode in (:none, :expand_energy_window) ||
            throw(ArgumentError("unknown domain adaptation policy"))
        expansions >= 0 || throw(ArgumentError("maximum_expansions cannot be negative"))
        nodes >= 3 || throw(ArgumentError("maximum_energy_nodes must be at least three"))
        isfinite(growth) && growth > 0 ||
            throw(ArgumentError("growth_fraction must be positive"))
        isfinite(tail) && 0 < tail < 1 ||
            throw(ArgumentError("tail_threshold must lie in (0,1)"))
        new(mode, Int(expansions), Int(nodes), Float64(growth), Float64(tail))
    end
end
