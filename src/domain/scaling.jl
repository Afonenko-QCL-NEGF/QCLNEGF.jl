"""
    ScaleSystem(; E₀=0.1u"eV", L₀=10u"nm", m_ref=0.067CODATA.m₀)

Construct the nondimensionalization used by every dense linear-algebra
operation.  The Poisson coupling is computed with Unitful before units are
stripped, so an inconsistent dimensional formula cannot silently pass.

See [Dimensionless scaling](@ref theory-scaling).
"""
function ScaleSystem(; E₀ = 0.1u"eV", L₀ = 10.0u"nm", m_ref = 0.067 * CODATA.m₀)
    E₀q = _energy(E₀)
    L₀q = _length(L₀)
    mrefq = _mass(m_ref)
    E₀q > 0u"eV" || throw(ArgumentError("E₀ must be positive"))
    L₀q > 0u"m" || throw(ArgumentError("L₀ must be positive"))
    mrefq > 0u"kg" || throw(ArgumentError("m_ref must be positive"))
    E₀_J = Float64(ustrip(u"J", uconvert(u"J", E₀q)))
    E₀_eV = _electronvolts(E₀q)
    L₀_m = _metres(L₀q)
    λ_P_q = CODATA.e^2 / (CODATA.ε₀ * uconvert(u"J", E₀q) * L₀q)
    λ_P = Float64(ustrip(Unitful.NoUnits, uconvert(Unitful.NoUnits, λ_P_q)))
    J₀ = uconvert(u"A/m^2", CODATA.e * E₀q / (CODATA.ħₑᵥ * L₀q^2))
    return ScaleSystem(E₀q, L₀q, mrefq, E₀_eV, E₀_J, L₀_m, λ_P, J₀)
end

"""
    scaled_value(value, scales, kind)

Convert a Unitful value to the dimensionless internal representation.
Supported `kind`s are `:energy`, `:length`, `:wavenumber`, `:field`,
`:sheet_density`, and `:volume_density`.

See [Unitful and the dimensionless core](@ref theory-scaling).
"""
function scaled_value(value, s::ScaleSystem, kind::Symbol)
    if kind === :energy
        return _electronvolts(value) / s.E₀_eV
    elseif kind === :length
        return _metres(value) / s.L₀_m
    elseif kind === :wavenumber
        return _inverse_metres(value) * s.L₀_m
    elseif kind === :field
        F₀ = s.E₀_eV / s.L₀_m
        return _volts_per_metre(value) / F₀
    elseif kind === :sheet_density
        return _per_square_metre(value) * s.L₀_m^2
    elseif kind === :volume_density
        return _per_cubic_metre(value) * s.L₀_m^3
    end
    throw(ArgumentError("unsupported scale kind: $kind"))
end

"""
Inverse of [`scaled_value`](@ref), returning a Unitful quantity.

See [Unitful and the dimensionless core](@ref theory-scaling) and the
[result-unit contract](@ref full-input-contract).
"""
function physical_value(value::Real, s::ScaleSystem, kind::Symbol)
    if kind === :energy
        return value * s.E₀
    elseif kind === :length
        return value * s.L₀
    elseif kind === :wavenumber
        return value / s.L₀
    elseif kind === :field
        return uconvert(u"V/m", value * s.E₀ / (CODATA.e * s.L₀))
    elseif kind === :sheet_density
        return uconvert(u"m^-2", value / s.L₀^2)
    elseif kind === :volume_density
        return uconvert(u"m^-3", value / s.L₀^3)
    end
    throw(ArgumentError("unsupported scale kind: $kind"))
end

"""
    current_density_A_per_cm2(value_A_per_m2)

Convert the solver's canonical current-density storage unit (A/m²) to the
single human-presentation unit used by reports, terminal, HTML, and runtime
views (A/cm²).

See [Production observability](https://github.com/AfonenkoA/QCLNEGFRunner.jl/blob/main/docs/src/user/results.md).
"""
function current_density_A_per_cm2(value_A_per_m2::Real)
    value = Float64(value_A_per_m2)
    isfinite(value) ||
        throw(ArgumentError("current density must be finite for presentation"))
    quantity = value * u"A/m^2"
    return Float64(ustrip(u"A/cm^2", uconvert(u"A/cm^2", quantity)))
end

_C̄(s::ScaleSystem) =
    Float64(ustrip(Unitful.NoUnits, uconvert(Unitful.NoUnits, C₀ / (s.E₀ * s.L₀^2))))

function _scaled_physics(p::PhysicalParameters, s::ScaleSystem)
    Lp = scaled_value(period_length(p), s, :length)
    F = scaled_value(p.F_bias, s, :field)
    Eᵖ = F * Lp
    return (;
        Lp,
        F,
        Eᵖ,
        Nᴰ²ᴰ = scaled_value(ionized_sheet_density(p), s, :sheet_density),
        ħωᴸᴼ = scaled_value(p.ħωᴸᴼ, s, :energy),
        q_s = scaled_value(p.q_s, s, :wavenumber),
        qᴸᴼ_s = scaled_value(p.qᴸᴼ_s, s, :wavenumber),
        z₀ = scaled_value(p.z₀, s, :length),
        E_ref = scaled_value(p.E_ref, s, :energy),
    )
end
