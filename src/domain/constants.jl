# Fundamental constants are stored with Unitful dimensions.  The numerical
# values follow the SI/CODATA set fixed by the mathematical specification.

const LengthQuantity = typeof(1.0u"m")
const EnergyQuantity = typeof(1.0u"eV")
const MassQuantity = typeof(1.0u"kg")
const TemperatureQuantity = typeof(1.0u"K")
const FieldQuantity = typeof(1.0u"V/m")
const WaveNumberQuantity = typeof(1.0u"m^-1")
const SheetDensityQuantity = typeof(1.0u"m^-2")
const VolumeDensityQuantity = typeof(1.0u"m^-3")
const MassDensityQuantity = typeof(1.0u"kg/m^3")
const SpeedQuantity = typeof(1.0u"m/s")
const VolumeQuantity = typeof(1.0u"m^3")

"""
Dimensioned fundamental constants used by all model constructors.

See [the complete input contract](@ref full-input-contract) and
[dimensionless scaling](@ref theory-scaling).
"""
struct FundamentalConstants
    e::typeof(1.0u"C")
    ħ::typeof(1.0u"J*s")
    ħₑᵥ::typeof(1.0u"eV*s")
    kᴮ::typeof(1.0u"J/K")
    kᴮₑᵥ::typeof(1.0u"eV/K")
    m₀::MassQuantity
    ε₀::typeof(1.0u"F/m")
    c::SpeedQuantity
end

"""
`CODATA` is the immutable SI/eV constants table used by every conversion:
exact `e` and `kᴮ`, the stated CODATA values of `ħ`, `m₀`, and `ε₀`.  The
eV copies are derived with Unitful from the SI fields and are not independent
rounded constants.

See [Fundamental constants](@ref full-input-contract) for units and the role
of every stored value.
"""
const CODATA = let
    electron_charge = 1.602176634e-19u"C"
    hbar_si = 1.054571817e-34u"J*s"
    boltzmann_si = 1.380649e-23u"J/K"
    FundamentalConstants(
        electron_charge,
        hbar_si,
        uconvert(u"eV*s", hbar_si),
        boltzmann_si,
        uconvert(u"eV/K", boltzmann_si),
        9.1093837139e-31u"kg",
        8.8541878128e-12u"F/m",
        299_792_458.0u"m/s",
    )
end

# Derive the kinetic prefactor from the one canonical constants table. The
# written specification may show a rounded decimal, but executable physics
# must not carry a second independently rounded value.
const C₀ = uconvert(u"eV*nm^2", CODATA.ħ^2 / (2 * CODATA.m₀))

_q(unit, value) = Float64(ustrip(unit, uconvert(unit, value))) * unit
_length(value) = _q(u"m", value)
_energy(value) = _q(u"eV", value)
_mass(value) = _q(u"kg", value)
_temperature(value) = _q(u"K", value)
_field(value) = _q(u"V/m", value)
_wavenumber(value) = _q(u"m^-1", value)
_sheetdensity(value) = _q(u"m^-2", value)
_volumedensity(value) = _q(u"m^-3", value)
_massdensity(value) = _q(u"kg/m^3", value)
_speed(value) = _q(u"m/s", value)
_volume(value) = _q(u"m^3", value)

_metres(value) = Float64(ustrip(u"m", value))
_electronvolts(value) = Float64(ustrip(u"eV", value))
_inverse_metres(value) = Float64(ustrip(u"m^-1", value))
_kelvin(value) = Float64(ustrip(u"K", value))
_volts_per_metre(value) = Float64(ustrip(u"V/m", value))
_per_square_metre(value) = Float64(ustrip(u"m^-2", value))
_per_cubic_metre(value) = Float64(ustrip(u"m^-3", value))
function _software_version()
    version = Base.pkgversion(@__MODULE__)
    return version === nothing ? "unknown" : string(version)
end
