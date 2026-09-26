"""Estimated one-particle spectral coverage before an interacting solve.

All energies are in the declared solver energy unit relative to E_ref. The
bound covers every eigenvalue at every represented radial momentum. It does
not claim to bound the interacting spectral tails; those are measured later.
"""
struct RepresentationCoverage
    discretization_id::String
    basis_id::String
    hartree_id::String
    energy_min::Float64
    energy_max::Float64
    h_min::Float64
    h_max::Float64
    shift_margin::Float64
    hilbert_margin::Float64
    required_min::Float64
    required_max::Float64
    minimum_by_k::Vector{Float64}
    maximum_by_k::Vector{Float64}
    covered::Bool
    source::Symbol
end

function _representation_digest(parts...)
    # IDs use numerical contents and their shape, never Julia's salted hash.
    io = IOBuffer()
    for part in parts
        print(io, size(part), ':')
        for x in part
            print(io, repr(x), ';')
        end
        print(io, '|')
    end
    return bytes2hex(sha256(take!(io)))
end

function representation_coverage(
    problem::NEGFProblem;
    hartree::AbstractVector{<:Real} = zeros(problem.numerical.N_z),
)
    all(isfinite, hartree) ||
        throw(ArgumentError("Hartree coverage estimate must be finite"))
    h = project_hamiltonians(problem, hartree)
    nk = length(problem.grids.κ)
    lows = zeros(nk)
    highs = zeros(nk)
    for m = 1:nk
        block = Matrix(view(h, m, :, :))
        norm(block-block') <= 1e-11*max(norm(block), 1.0) ||
            throw(ArgumentError("coverage Hamiltonian is not Hermitian"))
        lows[m], highs[m] = extrema(eigvals(Hermitian(block)))
    end
    sp = _scaled_physics(problem.physical, problem.scales)
    shift = abs(sp.Eᵖ) + (problem.scattering.LO ? abs(sp.ħωᴸᴼ) : 0.0)
    if problem.models.electron_electron.mode !== :none
        pole = sppa_parameters(problem.physical, problem.models.electron_electron)
        shift += _electronvolts(pole.pole_energy)/problem.scales.E₀_eV
    end
    hilbert = scaled_value(problem.numerical.M_E, problem.scales, :energy)
    lower, upper = extrema(problem.grids.ε)
    needed_low = minimum(lows)-shift-hilbert
    needed_high = maximum(highs)+shift+hilbert
    g = problem.grids
    return RepresentationCoverage(
        _representation_digest(g.x, g.wˣ, g.ε, g.wᴱ, g.κ, g.wᵏ, g.qᶻ, g.φ),
        _representation_digest(problem.basis.Φ, problem.basis.H₀, problem.basis.M⁻¹),
        _representation_digest(hartree),
        lower,
        upper,
        minimum(lows),
        maximum(highs),
        shift,
        hilbert,
        needed_low,
        needed_high,
        lows,
        highs,
        lower <= needed_low && upper >= needed_high,
        :hamiltonian_estimate,
    )
end

"""Propose a wider energy window while preserving (or refining) the old step.

This returns a new NumericalParameters value. The caller must rebuild grids,
basis-dependent kernels and caches and record a new discretization; no data
transfer or exact restart is implied. The node budget fails explicitly.
"""
function expand_energy_window(
    problem::NEGFProblem,
    coverage::RepresentationCoverage;
    maximum_nodes::Integer = 100_001,
)
    current = representation_coverage(problem)
    coverage.discretization_id == current.discretization_id ||
        throw(ArgumentError("coverage does not belong to the current discretization"))
    n = problem.numerical
    step = (coverage.energy_max-coverage.energy_min)/(n.N_E-1)
    lower = min(coverage.energy_min, coverage.required_min)
    upper = max(coverage.energy_max, coverage.required_max)
    nodes = ceil(Int, (upper-lower)/step)+1
    nodes <= maximum_nodes || throw(
        ArgumentError(
            "requested spectral coverage exceeds the declared energy-node budget",
        ),
    )
    changes = (
        E_min = (lower*problem.scales.E₀_eV+_electronvolts(problem.physical.E_ref))*u"eV",
        E_max = (upper*problem.scales.E₀_eV+_electronvolts(problem.physical.E_ref))*u"eV",
        N_E = nodes,
    )
    return NumericalParameters(;
        (
            name=>(
                hasproperty(changes, name) ? getproperty(changes, name) : getfield(n, name)
            ) for name in fieldnames(NumericalParameters)
        )...,
    )
end

"""Actual spectral/occupied tails and matrix spectral deficit of one state.

The spectral sum deficit is a matrix norm, never a lost-charge fraction.
No spectral matrix or particle population is renormalized by this check.
"""
function measured_representation_coverage(
    problem::NEGFProblem,
    scba::SCBAResult,
    hartree::AbstractVector{<:Real};
    energy_window_fraction = 0.05,
    momentum_window_fraction = 0.1,
)
    coverage = representation_coverage(problem; hartree)
    total_sc = _sum_selfenergies(scba.scattering, size(scba.green.Gᴿ))
    tails = _edge_metrics(
        problem,
        scba.green,
        total_sc;
        energy_window_fraction,
        momentum_window_fraction,
    )
    return (;
        estimate = coverage,
        source = :measured_state,
        tails,
        spectral_matrix_deficit = _spectral_sum_residual(problem, scba.green),
        energy_step = maximum(diff(problem.grids.ε)),
        minimum_resolved_broadening = _minimum_positive_broadening(total_sc),
    )
end

function _minimum_positive_broadening(family::SelfEnergyFamily)
    value = Inf
    for e in axes(family.Σᴿ, 1), m in axes(family.Σᴿ, 2)
        block = im*(_matrix_block(family.Σᵍ, e, m)-_matrix_block(family.Σˡ, e, m))
        for eigenvalue in eigvals(Hermitian((block+block')/2))
            eigenvalue > 0 && (value=min(value, eigenvalue))
        end
    end
    return isfinite(value) ? value : nothing
end

"""Exact finite-window weight of a normalized Lorentzian spectral line."""
function lorentzian_window_weight(lower::Real, upper::Real, centre::Real, width::Real)
    all(isfinite, (lower, upper, centre, width)) && upper > lower && width > 0 ||
        throw(ArgumentError("Lorentzian bounds and positive width must be finite"))
    return (atan((upper-centre)/width)-atan((lower-centre)/width))/π
end
