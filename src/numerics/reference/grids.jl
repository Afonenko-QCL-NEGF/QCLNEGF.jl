"""
    build_grids(physical, numerical, scales=ScaleSystem())

Build the cell-centred spatial grid, trapezoidal energy and longitudinal
momentum quadratures, annular radial momentum weights, and periodic angular
nodes.  Every stored value is dimensionless; the exact axis order is part of
the [array contract](@ref array-contracts).

The associated physical cutoffs and units are specified in
[the complete input contract](@ref full-input-contract).
"""
function build_grids(
    p::PhysicalParameters,
    n::NumericalParameters,
    s::ScaleSystem = ScaleSystem(),
)
    sp = _scaled_physics(p, s)
    Δx = sp.Lp / n.N_z
    x = [(g - 0.5) * Δx for g = 1:n.N_z]
    wˣ = fill(Δx, n.N_z)

    εmin = (_electronvolts(n.E_min) - _electronvolts(p.E_ref)) / s.E₀_eV
    εmax = (_electronvolts(n.E_max) - _electronvolts(p.E_ref)) / s.E₀_eV
    ε = collect(range(εmin, εmax; length = n.N_E))
    Δε = (εmax - εmin) / (n.N_E - 1)
    wᴱ = fill(Δε, n.N_E)
    wᴱ[[1, end]] ./= 2
    margin = scaled_value(n.M_E, s, :energy)
    trusted_energy = BitVector((ε .>= εmin + margin) .& (ε .<= εmax - margin))

    κmax = scaled_value(n.k_max, s, :wavenumber)
    κ = collect(range(0.0, κmax; length = n.N_k))
    boundaries = zeros(Float64, n.N_k + 1)
    boundaries[1] = 0.0
    for m = 2:n.N_k
        boundaries[m] = (κ[m-1] + κ[m]) / 2
    end
    boundaries[end] = κmax
    wᵏ = [(boundaries[m+1]^2 - boundaries[m]^2) / (4π) for m = 1:n.N_k]

    qmax = scaled_value(n.qz_max, s, :wavenumber)
    qᶻ = collect(range(-qmax, qmax; length = n.N_qz))
    Δq = 2qmax / (n.N_qz - 1)
    wᑫᶻ = fill(Δq, n.N_qz)
    wᑫᶻ[[1, end]] ./= 2
    φ = [2π * (r - 1) / n.N_φ for r = 1:n.N_φ]

    return ModelGrids(x, wˣ, ε, wᴱ, trusted_energy, κ, wᵏ, qᶻ, wᑫᶻ, φ)
end

"""
    build_profiles(physical, numerical, scales, grids)

Map the six epitaxial segments to the cell centres.  The donor profile is
renormalized on the actual grid so that `dot(wˣ,Nᴰ) == Nᴰ²ᴰ` to roundoff.

See [the reference design structure passport](@ref theory-reference2019) and
[the array contract](@ref array-contracts).
"""
function build_profiles(
    p::PhysicalParameters,
    n::NumericalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
)
    right_edges = cumsum(scaled_value(layer.d, s, :length) for layer in p.layers)
    N = n.N_z
    Eᶜ = zeros(Float64, N)
    mᶻᵣ = zeros(Float64, N)
    m_parallelᵣ = zeros(Float64, N)
    εᵣ = zeros(Float64, N)
    x_Al = zeros(Float64, N)
    layer_index = zeros(Int, N)
    χᴰ = zeros(Float64, N)
    Eref = _electronvolts(p.E_ref)
    for g = 1:N
        # Half-open layer convention: z_left <= z < z_right.  An exactly
        # coincident interface therefore belongs to the layer on its right.
        ℓ = min(searchsortedlast(right_edges, grids.x[g]) + 1, length(p.layers))
        layer = p.layers[ℓ]
        layer_index[g] = ℓ
        Eᶜ[g] = (_electronvolts(layer.Eᶜ) - Eref) / s.E₀_eV
        mᶻᵣ[g] = layer.mᶻᵣ
        m_parallelᵣ[g] = layer.m_parallelᵣ
        εᵣ[g] = layer.εᵣ
        x_Al[g] = layer.x_Al
        χᴰ[g] = layer.doped ? 1.0 : 0.0
    end
    normalization = dot(grids.wˣ, χᴰ)
    normalization > 0 || throw(ArgumentError("no doped grid cells were selected"))
    Nᴰ²ᴰ = scaled_value(ionized_sheet_density(p), s, :sheet_density)
    Nᴰ = Nᴰ²ᴰ .* χᴰ ./ normalization
    return MaterialProfiles(Eᶜ, mᶻᵣ, m_parallelᵣ, εᵣ, Nᴰ, x_Al, layer_index)
end
