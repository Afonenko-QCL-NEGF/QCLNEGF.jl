"""
Opt-in physical primitives with explicit dimensions and independent limiting tests.

SCBA integration is performed by numerical adapters, with an unchanged default
baseline. See `docs/src/physics/models.md` for each integration boundary.
"""
module PhysicalModelExtensions

using LinearAlgebra
using Unitful
import ..QCLDomain: CODATA

export AbstractScreeningModel,
    FixedScreening,
    DebyeScreening3D,
    ThomasFermiScreening2D,
    screening_wavenumber,
    screened_coulomb_potential,
    model_identity,
    KaneDispersion,
    kane_energy,
    kane_velocity,
    kane_dos_2d,
    kane_state_count_2d,
    kane_energy_bin_weights_2d,
    kane_kinetic_operator,
    HotLOPhononBath,
    LOPopulation,
    equilibrium_lo_population,
    lo_scattering_factors,
    steady_lo_population,
    advance_lo_population,
    lo_energy_balance,
    lo_population_temperature,
    BinaryCollisionChannel,
    ElectronCollisionModel,
    electron_collision_rhs,
    electron_collision_moments,
    optical_power_density,
    single_plasmon_pole

const Energy = typeof(1.0u"eV")
const WaveNumber = typeof(1.0u"m^-1")
const Rate = typeof(1.0u"s^-1")

# Construction tolerances for checking the conservation identities of supplied
# discrete Pauli channels. These are absolute comparison thresholds, not
# material energies, momentum cutoffs, broadenings, or collision rates. They
# can be overridden explicitly (including zero) and are recorded in identity.
const DEFAULT_BINARY_ENERGY_CONSERVATION_ATOL = 1e-12u"eV"
const DEFAULT_BINARY_MOMENTUM_CONSERVATION_ATOL = 1e-4u"m^-1"

function _physical(unit, value, label; positive = false)
    result = Float64(ustrip(unit, uconvert(unit, value)))
    valid = isfinite(result) && (positive ? result > 0 : result >= 0)
    valid || throw(
        ArgumentError("$label must be finite and $(positive ? "positive" : "nonnegative")"),
    )
    return result * unit
end

function _positive(value::Real, label)
    result = Float64(value)
    isfinite(result) && result > 0 ||
        throw(ArgumentError("$label must be finite and positive"))
    return result
end

"""A physical screening closure; the spatial dimension is part of its identity."""
abstract type AbstractScreeningModel end

"""Fixed inverse screening length. `dimension` must be 2 or 3, explicitly."""
struct FixedScreening <: AbstractScreeningModel
    wavenumber::WaveNumber
    relative_permittivity::Float64
    dimension::Int
    function FixedScreening(wavenumber, relative_permittivity::Real, dimension::Integer)
        dimension in (2, 3) || throw(ArgumentError("screening dimension must be 2 or 3"))
        new(
            _physical(u"m^-1", wavenumber, "screening wavenumber"),
            _positive(relative_permittivity, "relative permittivity"),
            Int(dimension),
        )
    end
end

"""Classical, homogeneous 3D Debye closure at explicitly specified density and T."""
struct DebyeScreening3D <: AbstractScreeningModel
    density::typeof(1.0u"m^-3")
    temperature::typeof(1.0u"K")
    relative_permittivity::Float64
    function DebyeScreening3D(density, temperature, relative_permittivity::Real)
        new(
            _physical(u"m^-3", density, "volume density"),
            _physical(u"K", temperature, "electron temperature"; positive = true),
            _positive(relative_permittivity, "relative permittivity"),
        )
    end
end

"""Finite-T long-wavelength 2D compressibility of one ideal parabolic band."""
struct ThomasFermiScreening2D <: AbstractScreeningModel
    density::typeof(1.0u"m^-2")
    temperature::typeof(1.0u"K")
    band_mass::typeof(1.0u"kg")
    degeneracy::Int
    relative_permittivity::Float64
    function ThomasFermiScreening2D(
        density,
        temperature,
        band_mass,
        degeneracy::Integer,
        relative_permittivity::Real,
    )
        degeneracy > 0 || throw(ArgumentError("degeneracy must be positive"))
        new(
            _physical(u"m^-2", density, "sheet density"),
            _physical(u"K", temperature, "electron temperature"; positive = true),
            _physical(u"kg", band_mass, "band mass"; positive = true),
            Int(degeneracy),
            _positive(relative_permittivity, "relative permittivity"),
        )
    end
end

screening_wavenumber(model::FixedScreening) = model.wavenumber
screening_wavenumber(model::DebyeScreening3D) = uconvert(
    u"m^-1",
    sqrt(
        CODATA.e^2 * model.density /
        (CODATA.ε₀ * model.relative_permittivity * CODATA.kᴮ * model.temperature),
    ),
)

function screening_wavenumber(model::ThomasFermiScreening2D)
    density_of_states = model.degeneracy * model.band_mass / (2π * CODATA.ħ^2)
    degeneracy_parameter = ustrip(
        Unitful.NoUnits,
        model.density / (density_of_states * CODATA.kᴮ * model.temperature),
    )
    f_at_band_edge = -expm1(-degeneracy_parameter)
    return uconvert(
        u"m^-1",
        CODATA.e^2 * density_of_states * f_at_band_edge /
        (2 * CODATA.ε₀ * model.relative_permittivity),
    )
end

_dimension(model::FixedScreening) = model.dimension
_dimension(::DebyeScreening3D) = 3
_dimension(::ThomasFermiScreening2D) = 2

"""
Continuum Fourier potential, with inverse measure d^d q/(2π)^d.
Returns J m² in 2D and J m³ in 3D. No subband form factor or normalization
area/volume is included. An unscreened q=0 singularity is rejected.
"""
function screened_coulomb_potential(model::AbstractScreeningModel, wavenumber)
    q = _physical(u"m^-1", wavenumber, "transfer wavenumber")
    κ = screening_wavenumber(model)
    q + κ > 0u"m^-1" ||
        throw(DomainError(q, "unscreened zero-momentum Coulomb singularity"))
    ε = CODATA.ε₀ * model.relative_permittivity
    if _dimension(model) == 2
        return uconvert(u"J*m^2", CODATA.e^2 / (2 * ε * (q + κ)))
    end
    return uconvert(u"J*m^3", CODATA.e^2 / (ε * (q^2 + κ^2)))
end

model_identity(model::FixedScreening) = (
    kind = :fixed_screening,
    revision = 1,
    dimension = model.dimension,
    wavenumber_m_inverse = ustrip(u"m^-1", model.wavenumber),
    relative_permittivity = model.relative_permittivity,
)
model_identity(model::DebyeScreening3D) = (
    kind = :debye_3d,
    revision = 1,
    density_m3 = ustrip(u"m^-3", model.density),
    temperature_K = ustrip(u"K", model.temperature),
    relative_permittivity = model.relative_permittivity,
)
model_identity(model::ThomasFermiScreening2D) = (
    kind = :thomas_fermi_2d,
    revision = 1,
    density_m2 = ustrip(u"m^-2", model.density),
    temperature_K = ustrip(u"K", model.temperature),
    mass_kg = ustrip(u"kg", model.band_mass),
    degeneracy = model.degeneracy,
    relative_permittivity = model.relative_permittivity,
)

"""
Winge et al. (2016), equations (8)–(9): static/long-wavelength single plasmon
pole with pair-continuum term (ℏ²q²/2m)². The frequency-dependent correlation
weight has units of energy and multiplies the bare quasi-2D Coulomb potential.
The supplied electron temperature is a plasmon closure, not lattice temperature.
"""
function single_plasmon_pole(
    model::Union{DebyeScreening3D,ThomasFermiScreening2D},
    transfer_wavenumber,
    band_mass,
    electron_temperature,
)
    q = _physical(
        u"m^-1",
        transfer_wavenumber,
        "representative transfer momentum";
        positive = true,
    )
    mass = _physical(u"kg", band_mass, "effective mass"; positive = true)
    temperature =
        _physical(u"K", electron_temperature, "electron temperature"; positive = true)
    temperature == model.temperature ||
        throw(ArgumentError("plasmon and screening temperatures must agree"))
    model isa ThomasFermiScreening2D &&
        mass != model.band_mass &&
        throw(ArgumentError("plasmon and screening masses must agree"))
    κ = screening_wavenumber(model)
    κ > 0u"m^-1" || throw(
        ArgumentError("the plasmon-pole approximation requires positive carrier density"),
    )
    ε = CODATA.ε₀ * model.relative_permittivity
    plasma_squared = if model isa ThomasFermiScreening2D
        CODATA.ħ^2 * CODATA.e^2 * model.density * q / (2ε * mass)
    else
        CODATA.ħ^2 * CODATA.e^2 * model.density / (ε * mass)
    end
    static_factor = model isa ThomasFermiScreening2D ? 1 + q / κ : 1 + q^2 / κ^2
    pair_energy = CODATA.ħ^2 * q^2 / (2mass)
    pole = uconvert(u"eV", sqrt(plasma_squared * static_factor + pair_energy^2))
    weight = uconvert(u"eV", plasma_squared / (2pole))
    occupation = inv(expm1(ustrip(Unitful.NoUnits, pole / (CODATA.kᴮₑᵥ * temperature))))
    isfinite(pole) &&
    pole > 0u"eV" &&
    isfinite(weight) &&
    weight > 0u"eV" &&
    isfinite(occupation) &&
    occupation >= 0 || throw(
        DomainError(
            (pole, weight, occupation),
            "plasmon closure exceeds finite numerical range",
        ),
    )
    return (
        plasma_energy = uconvert(u"eV", sqrt(plasma_squared)),
        pole_energy = pole,
        correlation_weight = weight,
        occupation = occupation,
        screening_wavenumber = κ,
        transfer_wavenumber = q,
        screening_dimension = _dimension(model),
    )
end

"""Isotropic, homogeneous Kane band: E(1 + αE) = ℏ²k²/(2m). α is supplied."""
struct KaneDispersion
    band_mass::typeof(1.0u"kg")
    nonparabolicity::typeof(1.0u"eV^-1")
    degeneracy::Int
    function KaneDispersion(band_mass, nonparabolicity, degeneracy::Integer)
        degeneracy > 0 || throw(ArgumentError("degeneracy must be positive"))
        new(
            _physical(u"kg", band_mass, "band mass"; positive = true),
            _physical(u"eV^-1", nonparabolicity, "nonparabolicity"),
            Int(degeneracy),
        )
    end
end

model_identity(model::KaneDispersion) = (
    kind = :homogeneous_kane,
    revision = 1,
    mass_kg = ustrip(u"kg", model.band_mass),
    nonparabolicity_eV_inverse = ustrip(u"eV^-1", model.nonparabolicity),
    degeneracy = model.degeneracy,
)

function _kane_from_parabolic(model::KaneDispersion, parabolic_energy)
    e = _physical(u"eV", parabolic_energy, "parabolic kinetic energy")
    x = ustrip(Unitful.NoUnits, model.nonparabolicity * e)
    # Rationalized root has the correct α→0 limit without cancellation.
    return 2e / (1 + sqrt(1 + 4x))
end

function kane_energy(model::KaneDispersion, wavenumber)
    k = _physical(u"m^-1", wavenumber, "wavevector magnitude")
    return _kane_from_parabolic(model, CODATA.ħ^2 * k^2 / (2 * model.band_mass))
end

"""Group speed (1/ℏ)dE/dk, consistent with `kane_energy`."""
function kane_velocity(model::KaneDispersion, wavenumber)
    k = _physical(u"m^-1", wavenumber, "wavevector magnitude")
    energy = kane_energy(model, k)
    return uconvert(
        u"m/s",
        CODATA.ħ * k / (model.band_mass * (1 + 2 * model.nonparabolicity * energy)),
    )
end

"""2D density of states per energy per area, including the declared degeneracy."""
function kane_dos_2d(model::KaneDispersion, kinetic_energy)
    e = _physical(u"eV", kinetic_energy, "kinetic energy")
    return uconvert(
        u"eV^-1*m^-2",
        model.degeneracy * model.band_mass / (2π * CODATA.ħ^2) *
        (1 + 2 * model.nonparabolicity * e),
    )
end

"""Exact number of states per area between the band minimum and E."""
function kane_state_count_2d(model::KaneDispersion, kinetic_energy)
    e = _physical(u"eV", kinetic_energy, "kinetic energy")
    return uconvert(
        u"m^-2",
        model.degeneracy * model.band_mass / (2π * CODATA.ħ^2) *
        e *
        (1 + model.nonparabolicity * e),
    )
end

"""Energy-cell measures ∫cell D(E)dE; momentum-space dk measures stay unchanged."""
function kane_energy_bin_weights_2d(model::KaneDispersion, edges)
    length(edges) >= 2 || throw(ArgumentError("at least two energy edges are required"))
    e = [_physical(u"eV", value, "energy edge") for value in edges]
    all(diff(e) .> 0u"eV") || throw(ArgumentError("energy edges must increase strictly"))
    return diff(kane_state_count_2d.(Ref(model), e))
end

"""
Apply the same dispersion to a positive Hermitian homogeneous kinetic operator
whose entries carry energy units. This is a spectral function of T, not a
prescription for a position-dependent mass or for transforming T + V together.
"""
function kane_kinetic_operator(model::KaneDispersion, kinetic::AbstractMatrix)
    size(kinetic, 1) == size(kinetic, 2) && !isempty(kinetic) ||
        throw(ArgumentError("kinetic operator must be nonempty and square"))
    matrix = ComplexF64.(ustrip.(u"eV", uconvert.(u"eV", kinetic)))
    all(isfinite, matrix) || throw(ArgumentError("kinetic operator must be finite"))
    ishermitian(matrix) || throw(ArgumentError("kinetic operator must be Hermitian"))
    eigenbasis = eigen(Hermitian(matrix))
    minimum(eigenbasis.values) >= 0 || throw(
        DomainError(
            minimum(eigenbasis.values),
            "kinetic operator must be positive semidefinite; no eigenvalue clipping is performed",
        ),
    )
    values = ustrip.(u"eV", _kane_from_parabolic.(Ref(model), eigenbasis.values .* u"eV"))
    return (eigenbasis.vectors * Diagonal(values) * eigenbasis.vectors') .* u"eV"
end

"""Dispersionless LO mode with explicit lattice temperature and decay time."""
struct HotLOPhononBath
    phonon_energy::Energy
    lattice_temperature::typeof(1.0u"K")
    decay_time::typeof(1.0u"s")
    function HotLOPhononBath(phonon_energy, lattice_temperature, decay_time)
        new(
            _physical(u"eV", phonon_energy, "LO energy"; positive = true),
            _physical(u"K", lattice_temperature, "lattice temperature"),
            _physical(u"s", decay_time, "LO decay time"; positive = true),
        )
    end
end

"""Occupation of one LO mode (or one explicitly declared uniform mode class)."""
struct LOPopulation
    occupation::Float64
    function LOPopulation(occupation::Real)
        isfinite(occupation) && occupation >= 0 ||
            throw(ArgumentError("LO occupation must be finite and nonnegative"))
        new(Float64(occupation))
    end
end

model_identity(model::HotLOPhononBath) = (
    kind = :lo_rate_balance,
    revision = 1,
    phonon_energy_eV = ustrip(u"eV", model.phonon_energy),
    lattice_temperature_K = ustrip(u"K", model.lattice_temperature),
    decay_time_s = ustrip(u"s", model.decay_time),
)

function equilibrium_lo_population(model::HotLOPhononBath)
    model.lattice_temperature == 0u"K" && return LOPopulation(0)
    x = ustrip(
        Unitful.NoUnits,
        model.phonon_energy / (CODATA.kᴮₑᵥ * model.lattice_temperature),
    )
    return LOPopulation(inv(expm1(x)))
end

lo_scattering_factors(population::LOPopulation) =
    (emission = population.occupation + 1, absorption = population.occupation)

"""Equivalent Bose temperature of a uniform LO population, distinct from lattice T."""
function lo_population_temperature(model::HotLOPhononBath, population::LOPopulation)
    population.occupation == 0 && return 0.0u"K"
    return uconvert(
        u"K",
        model.phonon_energy / (CODATA.kᴮₑᵥ * log1p(inv(population.occupation))),
    )
end

function _lo_rates(model::HotLOPhononBath, emission_rate, absorption_rate)
    emission = _physical(u"s^-1", emission_rate, "spontaneous emission coefficient")
    absorption = _physical(u"s^-1", absorption_rate, "absorption coefficient")
    decay = inv(model.decay_time)
    source = emission + equilibrium_lo_population(model).occupation * decay
    damping = absorption + decay - emission
    return (
        source = source,
        damping = damping,
        emission = emission,
        absorption = absorption,
    )
end

"""Stationary solution of dN/dt = a(N+1) - bN - (N-Nbath)/τ."""
function steady_lo_population(model::HotLOPhononBath, emission_rate, absorption_rate)
    rates = _lo_rates(model, emission_rate, absorption_rate)
    rates.damping > 0u"s^-1" || throw(
        DomainError(
            rates.damping,
            "no stable finite stationary LO population for frozen electronic rates",
        ),
    )
    return LOPopulation(ustrip(Unitful.NoUnits, rates.source / rates.damping))
end

"""Exact nonnegative step with frozen electronic rates; unstable growth is explicit."""
function advance_lo_population(
    model::HotLOPhononBath,
    population::LOPopulation,
    emission_rate,
    absorption_rate,
    timestep,
)
    dt = _physical(u"s", timestep, "timestep")
    rates = _lo_rates(model, emission_rate, absorption_rate)
    x = ustrip(Unitful.NoUnits, rates.damping * dt)
    source_dt = ustrip(Unitful.NoUnits, rates.source * dt)
    response = x == 0 ? 1.0 : -expm1(-x) / x
    return LOPopulation(population.occupation * exp(-x) + source_dt * response)
end

"""Per-mode electron→LO, LO→bath and stored-energy rates, all in watts."""
function lo_energy_balance(
    model::HotLOPhononBath,
    population::LOPopulation,
    emission_rate,
    absorption_rate,
)
    rates = _lo_rates(model, emission_rate, absorption_rate)
    electron_to_lo = uconvert(
        u"W",
        model.phonon_energy * (
            rates.emission * (population.occupation + 1) -
            rates.absorption * population.occupation
        ),
    )
    lo_to_bath = uconvert(
        u"W",
        model.phonon_energy *
        (population.occupation - equilibrium_lo_population(model).occupation) /
        model.decay_time,
    )
    return (
        electron_to_lo = electron_to_lo,
        lo_to_bath = lo_to_bath,
        stored_energy_rate = electron_to_lo - lo_to_bath,
    )
end

"""
One reversible 1+2 ↔ 3+4 fermion collision between distinct spin-resolved states.
The nonnegative coefficient in s⁻¹ must already include the microscopic matrix
element, normalization, and quadrature. The reverse reaction is implicit.
"""
struct BinaryCollisionChannel
    states::NTuple{4,Int}
    rate::Rate
    function BinaryCollisionChannel(states::NTuple{4,<:Integer}, rate)
        all(>(0), states) && length(unique(states)) == 4 || throw(
            ArgumentError("a binary channel needs four distinct positive state indices"),
        )
        new(Int.(states), _physical(u"s^-1", rate, "collision coefficient"))
    end
end

"""
On-shell Markov/Pauli population collision model on equal-weight discrete states.
Every channel conserves energy and both in-plane momentum components. No Hartree
or exchange energy shift is added. A provenance label for the supplied rates is
mandatory. Input arrays are copied into immutable tuples.
`energy_tolerance` and `momentum_tolerance` are absolute construction-check
tolerances, defaulting to `DEFAULT_BINARY_ENERGY_CONSERVATION_ATOL` and
`DEFAULT_BINARY_MOMENTUM_CONSERVATION_ATOL`; neither changes a supplied state.
"""
struct ElectronCollisionModel{N,C}
    energies::NTuple{N,Energy}
    wavevectors::NTuple{N,NTuple{2,WaveNumber}}
    channels::NTuple{C,BinaryCollisionChannel}
    rate_provenance::String
    energy_tolerance::Energy
    momentum_tolerance::WaveNumber
    function ElectronCollisionModel(
        energies,
        wavevectors,
        channels;
        rate_provenance::AbstractString,
        energy_tolerance = DEFAULT_BINARY_ENERGY_CONSERVATION_ATOL,
        momentum_tolerance = DEFAULT_BINARY_MOMENTUM_CONSERVATION_ATOL,
    )
        n = length(energies)
        n > 0 && length(wavevectors) == n || throw(ArgumentError("state axes must agree"))
        isempty(strip(rate_provenance)) &&
            throw(ArgumentError("rate provenance is required"))
        e = Tuple(
            Float64(ustrip(u"eV", uconvert(u"eV", value))) * u"eV" for value in energies
        )
        all(isfinite, e) || throw(ArgumentError("quasiparticle energies must be finite"))
        k = Tuple(_signed_wavevector(value) for value in wavevectors)
        c = Tuple(channels)
        all(channel -> channel isa BinaryCollisionChannel, c) ||
            throw(ArgumentError("all collision channels must be BinaryCollisionChannel"))
        etol = _physical(u"eV", energy_tolerance, "energy tolerance")
        ktol = _physical(u"m^-1", momentum_tolerance, "momentum tolerance")
        seen = Set{NTuple{4,Int}}()
        for channel in c
            all(index -> index <= n, channel.states) ||
                throw(ArgumentError("channel index exceeds state axis"))
            a, b, d, f = channel.states
            abs(e[a] + e[b] - e[d] - e[f]) <= etol ||
                throw(ArgumentError("binary collision does not conserve energy"))
            all(
                axis -> abs(k[a][axis] + k[b][axis] - k[d][axis] - k[f][axis]) <= ktol,
                1:2,
            ) ||
                throw(ArgumentError("binary collision does not conserve in-plane momentum"))
            pair1, pair2 = minmax(a, b), minmax(d, f)
            left, right = minmax(pair1, pair2)
            key = (left..., right...)
            key in seen && throw(
                ArgumentError(
                    "duplicate or reversed binary collision double-counts a channel",
                ),
            )
            push!(seen, key)
        end
        new{n,length(c)}(e, k, c, String(rate_provenance), etol, ktol)
    end
end

function _signed_wavevector(value)
    length(value) == 2 ||
        throw(ArgumentError("each in-plane wavevector needs two components"))
    result = Tuple(
        Float64(ustrip(u"m^-1", uconvert(u"m^-1", component))) * u"m^-1" for
        component in value
    )
    all(isfinite, result) || throw(ArgumentError("wavevector must be finite"))
    return result
end

model_identity(model::ElectronCollisionModel) = (
    kind = :on_shell_pauli_binary,
    revision = 1,
    energies_eV = ustrip.(u"eV", model.energies),
    wavevectors_m_inverse = map(k -> ustrip.(u"m^-1", k), model.wavevectors),
    channels = map(
        c -> (states = c.states, rate_s_inverse = ustrip(u"s^-1", c.rate)),
        model.channels,
    ),
    rate_provenance = model.rate_provenance,
    energy_tolerance_eV = ustrip(u"eV", model.energy_tolerance),
    momentum_tolerance_m_inverse = ustrip(u"m^-1", model.momentum_tolerance),
)

"""Collision-only time derivative of spin-resolved occupations, in s⁻¹."""
function electron_collision_rhs(model::ElectronCollisionModel, occupations)
    length(occupations) == length(model.energies) ||
        throw(ArgumentError("occupation axis mismatch"))
    f = Float64.(occupations)
    all(value -> isfinite(value) && 0 <= value <= 1, f) ||
        throw(ArgumentError("fermion occupations must be finite in [0,1]"))
    derivative = zeros(length(f)) .* u"s^-1"
    for channel in model.channels
        a, b, c, d = channel.states
        forward = f[a] * f[b] * (1 - f[c]) * (1 - f[d])
        backward = f[c] * f[d] * (1 - f[a]) * (1 - f[b])
        flux = channel.rate * (forward - backward)
        derivative[a] -= flux
        derivative[b] -= flux
        derivative[c] += flux
        derivative[d] += flux
    end
    return derivative
end

"""Independent moment projections of a supplied collision derivative."""
function electron_collision_moments(model::ElectronCollisionModel, derivative)
    length(derivative) == length(model.energies) ||
        throw(ArgumentError("derivative axis mismatch"))
    rates = uconvert.(u"s^-1", derivative)
    all(isfinite, rates) || throw(ArgumentError("collision derivative must be finite"))
    return (
        number_rate = sum(rates),
        energy_rate = uconvert(u"W", sum(model.energies .* rates)),
        momentum_rate = ntuple(
            axis -> uconvert(
                u"N",
                CODATA.ħ * sum(map(k -> k[axis], model.wavevectors) .* rates),
            ),
            2,
        ),
    )
end

"""
Cycle-averaged field→matter power density, ε₀ω Im(χ)|E₀|²/2, for the
e⁻ⁱωᵗ convention and peak (not RMS) field amplitude. Negative values mean
stimulated power delivered to the field. This diagnostic does not update the
carrier population and introduces neither saturation nor cavity feedback.
"""
function optical_power_density(susceptibility::Number, angular_frequency, peak_field)
    isfinite(susceptibility) || throw(ArgumentError("susceptibility must be finite"))
    ω = _physical(u"s^-1", angular_frequency, "angular frequency"; positive = true)
    field = _physical(u"V/m", peak_field, "peak field magnitude")
    return uconvert(u"W/m^3", CODATA.ε₀ * ω * imag(susceptibility) * field^2 / 2)
end

end # module PhysicalModelExtensions
