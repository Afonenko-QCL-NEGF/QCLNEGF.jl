"""
Explicit representative-momentum single-plasmon-pole approximation.

Enabled calculations require transfer momentum, electron temperature, effective
mass and carrier density. The carrier-density dimension follows the declared
screening dimension. No material or electronic temperature is guessed.
"""
struct ElectronElectronOptions
    mode::Symbol
    screening_dimension::Int
    transfer_wavenumber::Union{Nothing,WaveNumberQuantity}
    electron_temperature::Union{Nothing,TemperatureQuantity}
    effective_mass_ratio::Union{Nothing,Float64}
    carrier_density_si::Union{Nothing,Float64}
    include_exchange::Bool
    function ElectronElectronOptions(
        mode::Symbol,
        dimension::Integer,
        transfer_wavenumber,
        electron_temperature,
        effective_mass_ratio::Union{Nothing,Real},
        carrier_density_si::Union{Nothing,Real},
        include_exchange::Bool,
    )
        mode in (:none, :sppa_single_q) ||
            throw(ArgumentError("unknown electron-electron model"))
        dimension in (2, 3) || throw(ArgumentError("screening_dimension must be 2 or 3"))
        supplied = (
            transfer_wavenumber,
            electron_temperature,
            effective_mass_ratio,
            carrier_density_si,
        )
        if mode === :none
            all(isnothing, supplied) || throw(
                ArgumentError(
                    "disabled e-e model cannot hide supplied physical parameters",
                ),
            )
            return new(
                mode,
                Int(dimension),
                nothing,
                nothing,
                nothing,
                nothing,
                include_exchange,
            )
        end
        any(isnothing, supplied) &&
            throw(ArgumentError("enabled SPPA requires every physical parameter"))
        q, temperature =
            _wavenumber(transfer_wavenumber), _temperature(electron_temperature)
        mass, density = Float64(effective_mass_ratio), Float64(carrier_density_si)
        isfinite(q) && q > 0u"m^-1" ||
            throw(ArgumentError("transfer wavenumber must be finite and positive"))
        isfinite(temperature) && temperature > 0u"K" ||
            throw(ArgumentError("electron temperature must be finite and positive"))
        isfinite(mass) && mass > 0 ||
            throw(ArgumentError("effective mass ratio must be finite and positive"))
        isfinite(density) && density > 0 ||
            throw(ArgumentError("carrier density must be finite and positive"))
        new(mode, Int(dimension), q, temperature, mass, density, include_exchange)
    end
end

function ElectronElectronOptions(;
    mode::Symbol = :none,
    screening_dimension::Integer = 2,
    transfer_wavenumber = nothing,
    electron_temperature = nothing,
    effective_mass_ratio = nothing,
    carrier_density = nothing,
    include_exchange::Bool = true,
)
    supplied =
        (transfer_wavenumber, electron_temperature, effective_mass_ratio, carrier_density)
    if mode === :none
        all(isnothing, supplied) || throw(
            ArgumentError("disabled e-e model cannot hide supplied physical parameters"),
        )
        return ElectronElectronOptions(
            :none,
            screening_dimension,
            nothing,
            nothing,
            nothing,
            nothing,
            include_exchange,
        )
    end
    any(isnothing, supplied) && throw(
        ArgumentError(
            "sppa_single_q requires transfer_wavenumber, electron_temperature, effective_mass_ratio and carrier_density",
        ),
    )
    unit = screening_dimension == 2 ? u"m^-2" : u"m^-3"
    density = Float64(ustrip(unit, uconvert(unit, carrier_density)))
    return ElectronElectronOptions(
        mode,
        screening_dimension,
        transfer_wavenumber,
        electron_temperature,
        effective_mass_ratio,
        density,
        include_exchange,
    )
end

function electron_electron_identity(options::ElectronElectronOptions)
    options.mode === :none && return (mode = :none, revision = 1)
    return (
        mode = options.mode,
        revision = 1,
        screening_dimension = options.screening_dimension,
        transfer_wavenumber_m_inverse = _inverse_metres(options.transfer_wavenumber),
        electron_temperature_K = _kelvin(options.electron_temperature),
        effective_mass_ratio = options.effective_mass_ratio,
        carrier_density_si = options.carrier_density_si,
        include_exchange = options.include_exchange,
        approximation = :representative_momentum_equilibrium_plasmon,
    )
end
