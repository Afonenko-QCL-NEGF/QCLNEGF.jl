"""
    projected_spectrum(hamiltonian_eV, W, metric; solve_rows, observation_rows)

Diagonalize the Hamiltonian in the column space of `W`. The Hamiltonian is a
matrix of operator entries in eV, acting on coordinate samples. The positive
diagonal `metric` specifies the inner product of those samples and must make
the Hamiltonian self-adjoint. For quadrature-normalized amplitudes use ones.

With `solve_rows`, restrict both the Hamiltonian and the columns of `W` to
those coordinates before projecting and diagonalizing. With
`observation_rows`, select amplitude rows after solving, without renormalizing
or changing any energy. Observation rows must belong to the solve domain. Both selections are vectors or ranges of integer coordinates.

Returns ascending `energies_eV`, sampled `amplitudes`, their `metric`, and the
coordinate indices `rows`. The columns are orthonormal on the solve domain.
Their norm can decrease on a smaller observation domain. The decomposition uses
Float64 precision, rejecting weighted Hermiticity defects above `100eps(Float64)`
relative and singular values at or below `100eps(Float64) * σmax`.
"""
function projected_spectrum(
    hamiltonian_eV::AbstractMatrix{<:Number},
    W::AbstractMatrix{<:Number},
    metric::AbstractVector{<:Real};
    solve_rows::AbstractVector{<:Integer} = Base.OneTo(size(W, 1)),
    observation_rows::AbstractVector{<:Integer} = solve_rows,
)
    Base.require_one_based_indexing(hamiltonian_eV, W, metric)
    samples, rank = size(W)
    size(hamiltonian_eV) == (samples, samples) && length(metric) == samples ||
        throw(DimensionMismatch("H, W and the metric must share their coordinate axis"))
    0 < rank <= length(solve_rows) ||
        throw(ArgumentError("basis rank must fit the solve domain"))
    all(row -> !(row isa Bool) && row in axes(W, 1), solve_rows) &&
        allunique(solve_rows) ||
        throw(ArgumentError("solve rows must be unique coordinate indices"))
    !isempty(observation_rows) && allunique(observation_rows) ||
        throw(ArgumentError("observation rows must be nonempty and unique"))
    positions = Dict(row => index for (index, row) in enumerate(solve_rows))
    all(row -> !(row isa Bool) && haskey(positions, row), observation_rows) ||
        throw(ArgumentError("observation rows must belong to the solve domain"))
    all(w -> isfinite(w) && w > 0, metric) ||
        throw(ArgumentError("the coordinate metric must be finite and positive"))
    all(isfinite, hamiltonian_eV) && all(isfinite, W) ||
        throw(ArgumentError("H and W must be finite"))

    weights = sqrt.(Float64.(metric[solve_rows]))
    all(w -> isfinite(w) && w > 0, weights) ||
        throw(ArgumentError("the solve metric must be representable in Float64"))
    h_type = eltype(hamiltonian_eV) <: Real ? Float64 : ComplexF64
    w_type = eltype(W) <: Real ? Float64 : ComplexF64
    orthogonal_H = weights .* Matrix{h_type}(hamiltonian_eV[solve_rows, solve_rows]) ./
                   transpose(weights)
    all(isfinite, orthogonal_H) ||
        throw(ArgumentError("weighted Hamiltonian is nonfinite"))
    isapprox(orthogonal_H, orthogonal_H'; rtol = 100eps(Float64), atol = 0) ||
        throw(ArgumentError("H is not self-adjoint in the coordinate metric"))

    # SVD of M^(1/2)W provides an orthonormal basis without squaring its
    # condition number in a Gram matrix.
    weighted_basis = weights .* Matrix{w_type}(W[solve_rows, :])
    all(isfinite, weighted_basis) || throw(ArgumentError("weighted basis is nonfinite"))
    subspace = svd(weighted_basis)
    last(subspace.S) > 100eps(Float64) * first(subspace.S) ||
        throw(ArgumentError("restricted basis is rank deficient in the coordinate metric"))
    reduced_H = subspace.U' * orthogonal_H * subspace.U
    spectrum = eigen(Hermitian((reduced_H + reduced_H') / 2))
    amplitudes = (subspace.U * spectrum.vectors) ./ weights
    selection = [positions[row] for row in observation_rows]
    return (;
        energies_eV = spectrum.values,
        amplitudes = amplitudes[selection, :],
        metric = metric[observation_rows],
        rows = collect(observation_rows),
    )
end

"""
    localized_cell_spectrum(physical, numerical, scales, grids, profiles; localization)

Construct the localized basis on the declared multi-period window and solve
the projected central-cell Hamiltonian at k = 0 and zero Hartree potential.
The bias in `physical` applies after constructing the basis. Returns absolute
energies in eV, coordinates and quadrature weights in nm, and wavefunction
samples in nm⁻¹ᐟ². The basis and projected Hamiltonian are included for
structural diagnostics.
"""
function localized_cell_spectrum(
    physical::PhysicalParameters,
    numerical::NumericalParameters,
    scales::ScaleSystem,
    grids::ModelGrids,
    profiles::MaterialProfiles;
    localization::Symbol = :pzp,
)
    iszero(first(grids.κ)) ||
        throw(ArgumentError("the cell spectrum requires an exact k = 0 ring"))
    basis = build_basis(physical, numerical, scales, grids, profiles; localization)
    projected = project_hamiltonians(
        basis, grids, profiles, physical, scales, zeros(length(grids.x)),
    )
    spectrum = eigen(Hermitian(Matrix(view(projected, 1, :, :))))
    nm = scales.L₀_m * 1e9
    metric_nm = grids.wˣ .* nm
    amplitudes = (basis.Φ * spectrum.vectors) ./ sqrt.(metric_nm)
    energy_reference_eV = _electronvolts(physical.E_ref)
    return (;
        basis,
        hamiltonians_dimensionless = projected,
        energies_eV = spectrum.values .* scales.E₀_eV .+ energy_reference_eV,
        energy_reference_eV,
        z_nm = grids.x .* nm,
        metric_nm,
        amplitudes_per_sqrt_nm = amplitudes,
    )
end
