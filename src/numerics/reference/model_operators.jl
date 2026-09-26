# Roundoff allowance for a signed scalar event sum. Multiplies eps, the
# multiplication/summation count, and sum of absolute weighted terms; no
# fixed physical rate or floor is introduced.
const _LO_EVENT_ROUNDOFF_FACTOR = 64.0

"""Apply a declared fixed thermodynamic screening closure before kernel assembly."""
function resolve_physical_models(p::PhysicalParameters, models::PhysicalModelOptions)
    models.screening === :fixed && return p
    density = ionized_sheet_density(p)
    temperature = models.screening_temperature_K*u"K"
    mass = models.screening_mass_ratio*CODATA.m₀
    dos = p.g_s*mass/(2π*CODATA.ħ^2)
    filling = -expm1(-ustrip(Unitful.NoUnits, density/(dos*CODATA.kᴮ*temperature)))
    q2d = _wavenumber(CODATA.e^2*dos*filling/(2CODATA.ε₀*p.ε_s))
    q3d = _wavenumber(
        sqrt(CODATA.e^2*density/(period_length(p)*CODATA.ε₀*p.ε_s*CODATA.kᴮ*temperature)),
    )
    return PhysicalParameters(
        (
            name === :q_s ? q2d : name === :qᴸᴼ_s ? q3d : getfield(p, name) for
            name in fieldnames(PhysicalParameters)
        )...,
    )
end

"""Local transverse Kane root; longitudinal confinement and all dk weights remain explicit."""
@inline function transverse_kinetic_energy(parabolic::Real, alpha_scaled::Real)
    isfinite(parabolic) && parabolic >= 0 && isfinite(alpha_scaled) && alpha_scaled >= 0 ||
        throw(
            ArgumentError(
                "transverse kinetic energy and nonparabolicity must be finite and nonnegative",
            ),
        )
    return 2parabolic/(1+sqrt(1+4alpha_scaled*parabolic))
end

function _transverse_hamiltonian(problem::NEGFProblem, momentum::Real)
    models = problem.models
    if models.dispersion === :parabolic
        return _C̄(problem.scales)*momentum^2*problem.basis.M⁻¹
    end
    alpha = models.nonparabolicity_per_eV*problem.scales.E₀_eV
    local_energy = [
        transverse_kinetic_energy(_C̄(problem.scales)*momentum^2/m, alpha) for
        m in problem.profiles.m_parallelᵣ
    ]
    return problem.basis.Φ'*Diagonal(local_energy)*problem.basis.Φ
end

"""Uniform LO population and energy-balance coefficients for the same raw map."""
struct LOKineticState
    occupation::Float64
    emission_per_second::Float64
    absorption_per_second::Float64
    electron_power_W_per_m2::Float64
    bath_power_W_per_m2::Float64
    balance_residual::Float64
    closure::Symbol
end

function _thermal_lo_occupation(p::PhysicalParameters)
    x = _electronvolts(p.ħωᴸᴼ)/_electronvolts(CODATA.kᴮₑᵥ*p.Tᴸᴼ)
    return inv(expm1(x))
end

function _lo_outgoing_coefficient(problem::NEGFProblem, green::GreenState, outgoing)
    rate = 0.0
    absolute = 0.0
    for e in axes(green.Gᴿ, 1), m in axes(green.Gᴿ, 2)
        block = real(tr(_matrix_block(outgoing, e, m)*_matrix_block(green.Gˡ, e, m)))
        weight = problem.grids.wᴱ[e]*problem.grids.wᵏ[m]
        rate += weight*block
        absolute += weight*abs(block)
    end
    operation_count =
        length(problem.grids.ε)*length(problem.grids.κ)*problem.numerical.N_b^2
    allowance = _LO_EVENT_ROUNDOFF_FACTOR*eps(Float64)*operation_count*absolute
    rate >= -allowance || throw(
        DomainError(
            rate,
            "negative electron-phonon transition rate; LO kinetic closure unusable",
        ),
    )
    # Only the declared scalar roundoff allowance is applied, never PSD clipping.
    rate = max(0.0, rate)
    flux_scale = problem.scales.E₀_eV/(ustrip(u"eV*s", CODATA.ħₑᵥ)*problem.scales.L₀_m^2)
    return problem.physical.g_s*rate/(2π)*flux_scale
end

function lo_kinetic_state(
    problem::NEGFProblem,
    green::GreenState,
    emission_outgoing,
    absorption_outgoing,
)
    models = problem.models
    emission_area = _lo_outgoing_coefficient(problem, green, emission_outgoing)
    absorption_area = _lo_outgoing_coefficient(problem, green, absorption_outgoing)
    mode_sheet = models.lo_mode_density_per_m3*_metres(period_length(problem.physical))
    bath = _thermal_lo_occupation(problem.physical)
    if models.lo_population === :rate_balance
        emission = emission_area/mode_sheet
        absorption = absorption_area/mode_sheet
        decay = inv(models.lo_decay_ps*1e-12)
        damping = absorption+decay-emission
        damping > 0 || throw(
            DomainError(damping, "uniform hot LO closure has no stable finite population"),
        )
        population = (emission+bath*decay)/damping
        bath_power =
            _electronvolts(problem.physical.ħωᴸᴼ) *
            _electron_charge() *
            mode_sheet *
            (population-bath) *
            decay
    else
        population =
            models.lo_population === :fixed_occupation ? models.lo_fixed_occupation : bath
        emission = NaN
        absorption = NaN
        bath_power = NaN
    end
    power =
        _electronvolts(problem.physical.ħωᴸᴼ) *
        _electron_charge() *
        (emission_area*(population+1)-absorption_area*population)
    balance =
        isfinite(bath_power) ?
        abs(power-bath_power)/max(abs(power), abs(bath_power), 1e-30) : NaN
    return LOKineticState(
        population,
        emission,
        absorption,
        power,
        bath_power,
        balance,
        models.lo_population,
    )
end

function _lo_population(problem::NEGFProblem, green, emission_outgoing, absorption_outgoing)
    model = problem.models
    model.lo_population === :thermal && return _thermal_lo_occupation(problem.physical)
    model.lo_population === :fixed_occupation && return model.lo_fixed_occupation
    return lo_kinetic_state(problem, green, emission_outgoing, absorption_outgoing).occupation
end
