"""Uniform validation failure at the physical/numerical scattering boundary."""
struct ScatteringValidationError <: Exception
    mechanism::Symbol
    field::String
    message::String
end

function Base.showerror(io::IO, error::ScatteringValidationError)
    print(
        io,
        "scattering validation error [",
        error.mechanism,
        "] at ",
        error.field,
        ": ",
        error.message,
    )
end

_scattering_error(mechanism, field, message) =
    throw(ScatteringValidationError(Symbol(mechanism), String(field), String(message)))

const _SCATTERING_MECHANISM_ORDER = (:LO, :acoustic, :impurity, :IFR, :alloy)

"""
Physical-model abstraction for one microscopic scattering mechanism.
Concrete models contain physical quantities only; grid choices, interpolation,
contraction order, threading, and storage do not belong here.
"""
abstract type AbstractScatteringPhysicalModel end

"""Bulk-like screened Fröhlich LO-phonon physical parameters."""
struct LOPhononModel <: AbstractScatteringPhysicalModel
    phonon_energy::EnergyQuantity
    temperature::TemperatureQuantity
    static_relative_permittivity::Float64
    high_frequency_relative_permittivity::Float64
    screening_wavenumber::WaveNumberQuantity
end


"""Elastic acoustic deformation-potential physical parameters."""
struct AcousticPhononModel <: AbstractScatteringPhysicalModel
    deformation_potential::EnergyQuantity
    temperature::TemperatureQuantity
    mass_density::MassDensityQuantity
    sound_velocity::SpeedQuantity
end

"""Screened ionized-donor physical parameters."""
struct IonizedImpurityModel <: AbstractScatteringPhysicalModel
    donor_sheet_density::SheetDensityQuantity
    ionization_fraction::Float64
    static_relative_permittivity::Float64
    screening_wavenumber::WaveNumberQuantity
end

"""Gaussian independent-interface roughness physical parameters."""
struct InterfaceRoughnessModel <: AbstractScatteringPhysicalModel
    rms_height::LengthQuantity
    correlation_length::LengthQuantity
    interfaces::Vector{LengthQuantity}
end

"""Short-range random-alloy physical parameters."""
struct AlloyDisorderModel <: AbstractScatteringPhysicalModel
    alloy_potential::EnergyQuantity
    primitive_cell_volume::VolumeQuantity
end

function _scattering_quantity(mechanism, field, converter, value)
    try
        return converter(value)
    catch error
        _scattering_error(
            mechanism,
            field,
            "invalid physical quantity: $(sprint(showerror, error))",
        )
    end
end

function LOPhononModel(;
    phonon_energy,
    temperature,
    static_relative_permittivity::Real,
    high_frequency_relative_permittivity::Real,
    screening_wavenumber,
)
    model = LOPhononModel(
        _scattering_quantity(:LO, "phonon_energy", _energy, phonon_energy),
        _scattering_quantity(:LO, "temperature", _temperature, temperature),
        Float64(static_relative_permittivity),
        Float64(high_frequency_relative_permittivity),
        _scattering_quantity(
            :LO,
            "screening_wavenumber",
            _wavenumber,
            screening_wavenumber,
        ),
    )
    return validate_scattering_model(model)
end

function AcousticPhononModel(;
    deformation_potential,
    temperature,
    mass_density,
    sound_velocity,
)
    model = AcousticPhononModel(
        _scattering_quantity(
            :acoustic,
            "deformation_potential",
            _energy,
            deformation_potential,
        ),
        _scattering_quantity(:acoustic, "temperature", _temperature, temperature),
        _scattering_quantity(:acoustic, "mass_density", _massdensity, mass_density),
        _scattering_quantity(:acoustic, "sound_velocity", _speed, sound_velocity),
    )
    return validate_scattering_model(model)
end

function IonizedImpurityModel(;
    donor_sheet_density,
    ionization_fraction::Real,
    static_relative_permittivity::Real,
    screening_wavenumber,
)
    model = IonizedImpurityModel(
        _scattering_quantity(
            :impurity,
            "donor_sheet_density",
            _sheetdensity,
            donor_sheet_density,
        ),
        Float64(ionization_fraction),
        Float64(static_relative_permittivity),
        _scattering_quantity(
            :impurity,
            "screening_wavenumber",
            _wavenumber,
            screening_wavenumber,
        ),
    )
    return validate_scattering_model(model)
end


function InterfaceRoughnessModel(; rms_height, correlation_length, interfaces)
    converted_interfaces = LengthQuantity[
        _scattering_quantity(:IFR, "interfaces[$index]", _length, value) for
        (index, value) in pairs(interfaces)
    ]
    model = InterfaceRoughnessModel(
        _scattering_quantity(:IFR, "rms_height", _length, rms_height),
        _scattering_quantity(:IFR, "correlation_length", _length, correlation_length),
        converted_interfaces,
    )
    return validate_scattering_model(model)
end

function AlloyDisorderModel(; alloy_potential, primitive_cell_volume)
    model = AlloyDisorderModel(
        _scattering_quantity(:alloy, "alloy_potential", _energy, alloy_potential),
        _scattering_quantity(
            :alloy,
            "primitive_cell_volume",
            _volume,
            primitive_cell_volume,
        ),
    )
    return validate_scattering_model(model)
end

"""Stable configuration/provenance identifier for a scattering model."""
scattering_id(::LOPhononModel) = :LO
scattering_id(::AcousticPhononModel) = :acoustic
scattering_id(::IonizedImpurityModel) = :impurity
scattering_id(::InterfaceRoughnessModel) = :IFR
scattering_id(::AlloyDisorderModel) = :alloy

"""
Any change to the enabled physical mechanism set is E3. This is distinct from
the E0--E2 numerical implementation selected for a fixed mechanism set.
"""
scattering_model_role(::AbstractScatteringPhysicalModel) = :physical_model
"""Return `:E0` for identical physical models and `:E3` for any change."""
function scattering_selection_evidence_class(left, right)
    allunique(scattering_id.(left)) || return :E3
    allunique(scattering_id.(right)) || return :E3
    right_by_id = Dict(scattering_id(model) => model for model in right)
    length(left) == length(right) && all(
        model ->
            haskey(right_by_id, scattering_id(model)) &&
            _same_scattering_model(model, right_by_id[scattering_id(model)]),
        left,
    ) ? :E0 : :E3
end

function _validate_positive_quantity(mechanism, field, value, zero)
    value > zero || _scattering_error(mechanism, field, "must be positive")
    return value
end

"""Fail-closed validation shared by every scattering physical model."""
function validate_scattering_model(model::LOPhononModel)
    mechanism = scattering_id(model)
    _validate_positive_quantity(mechanism, "phonon_energy", model.phonon_energy, 0u"eV")
    _validate_positive_quantity(mechanism, "temperature", model.temperature, 0u"K")
    model.static_relative_permittivity > 0 ||
        _scattering_error(mechanism, "static_relative_permittivity", "must be positive")
    0 < model.high_frequency_relative_permittivity < model.static_relative_permittivity ||
        _scattering_error(
            mechanism,
            "high_frequency_relative_permittivity",
            "must be positive and smaller than static permittivity",
        )
    _validate_positive_quantity(
        mechanism,
        "screening_wavenumber",
        model.screening_wavenumber,
        0u"m^-1",
    )
    return model
end

function validate_scattering_model(model::AcousticPhononModel)
    mechanism = scattering_id(model)
    _validate_positive_quantity(
        mechanism,
        "deformation_potential",
        model.deformation_potential,
        0u"eV",
    )
    _validate_positive_quantity(mechanism, "temperature", model.temperature, 0u"K")
    _validate_positive_quantity(mechanism, "mass_density", model.mass_density, 0u"kg/m^3")
    _validate_positive_quantity(mechanism, "sound_velocity", model.sound_velocity, 0u"m/s")
    return model
end

function validate_scattering_model(model::IonizedImpurityModel)
    mechanism = scattering_id(model)
    _validate_positive_quantity(
        mechanism,
        "donor_sheet_density",
        model.donor_sheet_density,
        0u"m^-2",
    )
    0 < model.ionization_fraction <= 1 ||
        _scattering_error(mechanism, "ionization_fraction", "must lie in (0,1]")
    model.static_relative_permittivity > 0 ||
        _scattering_error(mechanism, "static_relative_permittivity", "must be positive")
    _validate_positive_quantity(
        mechanism,
        "screening_wavenumber",
        model.screening_wavenumber,
        0u"m^-1",
    )
    return model
end

function validate_scattering_model(model::InterfaceRoughnessModel)
    mechanism = scattering_id(model)
    model.rms_height > 0u"m" || _scattering_error(
        mechanism,
        "rms_height",
        "must be positive for an enabled interface-roughness mechanism",
    )
    _validate_positive_quantity(
        mechanism,
        "correlation_length",
        model.correlation_length,
        0u"m",
    )
    isempty(model.interfaces) && _scattering_error(
        mechanism,
        "interfaces",
        "at least one interface is required for an enabled mechanism",
    )
    issorted(model.interfaces) ||
        _scattering_error(mechanism, "interfaces", "must be sorted")
    allunique(model.interfaces) ||
        _scattering_error(mechanism, "interfaces", "must not contain duplicates")
    return model
end

function validate_scattering_model(model::AlloyDisorderModel)
    mechanism = scattering_id(model)
    _validate_positive_quantity(mechanism, "alloy_potential", model.alloy_potential, 0u"eV")
    _validate_positive_quantity(
        mechanism,
        "primitive_cell_volume",
        model.primitive_cell_volume,
        0u"m^3",
    )
    return model
end

function validate_scattering_model(
    model::AbstractScatteringPhysicalModel,
    physical::PhysicalParameters,
)
    validate_scattering_model(model)
    return model
end

function validate_scattering_model(
    model::IonizedImpurityModel,
    physical::PhysicalParameters,
)
    validate_scattering_model(model)
    any(layer -> layer.doped, physical.layers) || _scattering_error(
        scattering_id(model),
        "doping_profile",
        "at least one epitaxial layer must carry donors",
    )
    return model
end

function validate_scattering_model(model::AlloyDisorderModel, physical::PhysicalParameters)
    validate_scattering_model(model)
    any(layer -> 0 < layer.x_Al < 1, physical.layers) || _scattering_error(
        scattering_id(model),
        "alloy_profile",
        "at least one layer must have an alloy fraction strictly inside (0,1)",
    )
    return model
end

function _lo_model(physical::PhysicalParameters)
    LOPhononModel(physical.ħωᴸᴼ, physical.Tᴸᴼ, physical.ε_s, physical.ε_∞, physical.qᴸᴼ_s)
end

function _acoustic_model(physical::PhysicalParameters)
    AcousticPhononModel(physical.Ξ, physical.Tᴸ, physical.ρ_m, physical.v_s)
end

function _impurity_model(physical::PhysicalParameters)
    IonizedImpurityModel(physical.N_dop²ᴰ, physical.f_ion, physical.ε_s, physical.q_s)
end

function _ifr_model(physical::PhysicalParameters)
    InterfaceRoughnessModel(physical.Δᴵᶠᴿ, physical.Λᴵᶠᴿ, copy(physical.interfaces))
end

function _alloy_model(physical::PhysicalParameters)
    physical.ΔV_alloy === nothing && _scattering_error(
        :alloy,
        "alloy_potential",
        "is required when alloy scattering is enabled",
    )
    physical.Ω₀ === nothing && _scattering_error(
        :alloy,
        "primitive_cell_volume",
        "is required when alloy scattering is enabled",
    )
    return AlloyDisorderModel(physical.ΔV_alloy, physical.Ω₀)
end

"""Build and uniformly validate the physical mechanism set selected by the
public `ScatteringOptions` value object."""
function scattering_models(physical::PhysicalParameters, options::ScatteringOptions)
    models = AbstractScatteringPhysicalModel[]
    options.LO && push!(models, _lo_model(physical))
    options.acoustic && push!(models, _acoustic_model(physical))
    options.impurity && push!(models, _impurity_model(physical))
    options.IFR && push!(models, _ifr_model(physical))
    options.alloy && push!(models, _alloy_model(physical))
    foreach(model -> validate_scattering_model(model, physical), models)
    return models
end

"""Project a validated physical-model collection to its mechanism selection.

`ScatteringOptions` is the stable selection value object shared by the direct
API and resolved run configuration. Physical parameters remain in the typed
models and are never reconstructed by this projection.
"""
function scattering_options(models::AbstractVector{<:AbstractScatteringPhysicalModel})
    identifiers = scattering_id.(models)
    allunique(identifiers) ||
        _scattering_error(:set, "mechanisms", "mechanism identifiers must be unique")
    return ScatteringOptions(
        LO = :LO in identifiers,
        acoustic = :acoustic in identifiers,
        impurity = :impurity in identifiers,
        IFR = :IFR in identifiers,
        alloy = :alloy in identifiers,
    )
end

_model_from_physical(::LOPhononModel, physical) = _lo_model(physical)
_model_from_physical(::AcousticPhononModel, physical) = _acoustic_model(physical)
_model_from_physical(::IonizedImpurityModel, physical) = _impurity_model(physical)
_model_from_physical(::InterfaceRoughnessModel, physical) = _ifr_model(physical)
_model_from_physical(::AlloyDisorderModel, physical) = _alloy_model(physical)
_model_from_physical(model::AbstractScatteringPhysicalModel, physical) = _scattering_error(
    scattering_id(model),
    "physical_parameters",
    "no canonical parameter projection is registered for $(typeof(model))",
)

function _models_match_physical(models, physical)
    return all(
        model -> _same_scattering_model(model, _model_from_physical(model, physical)),
        models,
    )
end

_same_scattering_model(left::LOPhononModel, right::LOPhononModel) =
    left.phonon_energy == right.phonon_energy &&
    left.temperature == right.temperature &&
    left.static_relative_permittivity == right.static_relative_permittivity &&
    left.high_frequency_relative_permittivity ==
    right.high_frequency_relative_permittivity &&
    left.screening_wavenumber == right.screening_wavenumber

_same_scattering_model(left::AcousticPhononModel, right::AcousticPhononModel) =
    left.deformation_potential == right.deformation_potential &&
    left.temperature == right.temperature &&
    left.mass_density == right.mass_density &&
    left.sound_velocity == right.sound_velocity

_same_scattering_model(left::IonizedImpurityModel, right::IonizedImpurityModel) =
    left.donor_sheet_density == right.donor_sheet_density &&
    left.ionization_fraction == right.ionization_fraction &&
    left.static_relative_permittivity == right.static_relative_permittivity &&
    left.screening_wavenumber == right.screening_wavenumber

_same_scattering_model(left::InterfaceRoughnessModel, right::InterfaceRoughnessModel) =
    left.rms_height == right.rms_height &&
    left.correlation_length == right.correlation_length &&
    left.interfaces == right.interfaces

_same_scattering_model(left::AlloyDisorderModel, right::AlloyDisorderModel) =
    left.alloy_potential == right.alloy_potential &&
    left.primitive_cell_volume == right.primitive_cell_volume

_same_scattering_model(
    ::AbstractScatteringPhysicalModel,
    ::AbstractScatteringPhysicalModel,
) = false
