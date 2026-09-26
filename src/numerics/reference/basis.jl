function _phase_fix!(Φ::AbstractMatrix{<:Complex})
    for a in axes(Φ, 2)
        g = argmax(abs.(view(Φ, :, a)))
        phase = angle(Φ[g, a])
        Φ[:, a] .*= exp(-im * phase)
        real(Φ[g, a]) < 0 && (Φ[:, a] .*= -1)
    end
    return Φ
end

function _repeat_profile(v::AbstractVector, count::Int)
    return repeat(v, count)
end

function _build_unlocalized_basis(
    p::PhysicalParameters,
    n::NumericalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    profiles::MaterialProfiles,
    Hwin::AbstractMatrix,
    Egrid₀::Real,
    Δx::Real,
)
    rows = (n.P_basis*n.N_z+1):((n.P_basis+1)*n.N_z)
    Hcell = Matrix(Hwin[rows, rows])
    cell_eigen = eigen(Hermitian(Hcell / Egrid₀))
    energies = cell_eigen.values[1:n.N_b] .* Egrid₀
    B = ComplexF64.(cell_eigen.vectors[:, 1:n.N_b])
    r_eigen = 0.0
    Hnorm = opnorm(Hcell)
    for j = 1:n.N_b
        bj = view(B, :, j)
        Ej = energies[j]
        r_eigen = max(
            r_eigen,
            norm(Hcell * bj - Ej * bj) / (Hnorm * norm(bj) + abs(Ej) * norm(bj)),
        )
    end

    S = Hermitian((B' * B + (B' * B)') / 2)
    eigS = eigen(S)
    smin, smax = extrema(eigS.values)
    κS = smax / smin
    smin > 0 || throw(ArgumentError("unlocalized cell eigenbasis is rank deficient"))
    Sinvhalf = eigS.vectors * Diagonal(1 ./ sqrt.(eigS.values)) * eigS.vectors'
    Φ = ComplexF64.(B * Sinvhalf)
    _phase_fix!(Φ)
    r_orth = norm(Φ' * Φ - I)

    H₀ = Matrix(Φ' * Hcell * Φ)
    Z = Matrix(Φ' * Diagonal(grids.x) * Φ)
    M⁻¹ = Matrix(Φ' * Diagonal(1 ./ profiles.m_parallelᵣ) * Φ)
    right_rows = (((n.P_basis+1)*n.N_z+1):((n.P_basis+2)*n.N_z))
    Hright = Matrix(Hwin[rows, right_rows])
    T₊ = Matrix(Φ' * Hright * Φ)
    T₋ = Matrix(T₊')
    χ = Φ ./ sqrt(Δx)
    centres = real.(diag(Z))
    second = real.(diag(Φ' * Diagonal(grids.x .^ 2) * Φ))
    spreads = sqrt.(max.(second .- centres .^ 2, 0.0))
    return BasisData(
        :none,
        Φ,
        χ,
        centres,
        spreads,
        collect(eigS.values),
        κS,
        H₀,
        Z,
        M⁻¹,
        T₊,
        T₋,
        collect(energies),
        Float64(Egrid₀),
        r_orth,
        r_eigen,
        0.0,
        0.0,
        0.0,
    )
end

"""
    build_localized_basis(physical, numerical, scales, grids, profiles)

Construct the fixed seven-period (or configured odd-window) zero-field
eigenspace, diagonalize `PzP`, select the central-period columns, perform the
cell-local Löwdin orthogonalization, and project `H₀`, `Z`, `M⁻¹`, and `T±`.

All returned matrices are dense.  `Φ` has shape `(N_z,N_b)` and is orthonormal
in the discrete cell metric; `χ=Φ/sqrt(Δx)` represents the continuous
dimensionless envelopes.  The function stops rather than pseudo-inverting a
rank-deficient cell overlap.

Implements [EQ-BASIS-001](@ref eq-pzp-lowdin).
See [BenDaniel–Duke and the fixed PzP basis](@ref theory-basis) for the
projection, localization, and Löwdin steps.
"""
function build_localized_basis(
    p::PhysicalParameters,
    n::NumericalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    profiles::MaterialProfiles;
    localization::Symbol = :pzp,
)
    localization in (:pzp, :none, :real_space) ||
        throw(ArgumentError("unsupported cell-basis localization $localization"))
    Np = 2n.P_basis + 1
    Nsub = Np * n.N_b
    Nwin = Np * n.N_z
    Nsub <= Nwin ||
        throw(ArgumentError("basis subspace cannot exceed the coordinate window"))
    Δx = grids.wˣ[1]
    Ewin = _repeat_profile(profiles.Eᶜ, Np)
    mwin = _repeat_profile(profiles.mᶻᵣ, Np)
    Hwin = build_bdd_hamiltonian(Ewin, mwin, Δx, _C̄(s))

    # Exact scalar rescaling from the coordinate diagonal scale.  It is close
    # to max_g(t_{g-1/2}+t_{g+1/2}) and does not double-count off-diagonals.
    Egrid₀ = maximum(abs, diag(Hwin))
    if localization === :real_space
        n.N_b == n.N_z || throw(
            ArgumentError("real_space requires N_b == N_z; no grid node may be discarded"),
        )
        rows = (n.P_basis*n.N_z+1):((n.P_basis+1)*n.N_z)
        right_rows = ((n.P_basis+1)*n.N_z+1):((n.P_basis+2)*n.N_z)
        Φ = Matrix{ComplexF64}(I, n.N_z, n.N_z)
        H₀ = ComplexF64.(Hwin[rows, rows])
        T₊ = ComplexF64.(Hwin[rows, right_rows])
        return BasisData(
            :real_space,
            Φ,
            Φ ./ sqrt(Δx),
            copy(grids.x),
            zeros(n.N_z),
            ones(n.N_z),
            1.0,
            H₀,
            Matrix{ComplexF64}(Diagonal(grids.x)),
            Matrix{ComplexF64}(Diagonal(1 ./ profiles.m_parallelᵣ)),
            T₊,
            Matrix(T₊'),
            eigvals(Hermitian(H₀)),
            Float64(Egrid₀),
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        )
    end
    localization === :none &&
        return _build_unlocalized_basis(p, n, s, grids, profiles, Hwin, Egrid₀, Δx)
    eigenwindow = eigen(Hermitian(Hwin / Egrid₀))
    eigenvalues_window = eigenwindow.values[1:Nsub] .* Egrid₀
    C = ComplexF64.(eigenwindow.vectors[:, 1:Nsub])
    r_eigen = 0.0
    Hnorm = opnorm(Hwin)
    for j = 1:Nsub
        cj = view(C, :, j)
        Ej = eigenvalues_window[j]
        r_eigen = max(
            r_eigen,
            norm(Hwin * cj - Ej * cj) / (Hnorm * norm(cj) + abs(Ej) * norm(cj)),
        )
    end

    xwin = [(G - 0.5) * Δx for G = 1:Nwin]
    Zsub = C' * (reshape(xwin, :, 1) .* C)
    zdiag = eigen(Hermitian((Zsub + Zsub') / 2))
    order = sortperm(zdiag.values)
    centres_all = zdiag.values[order]
    U = zdiag.vectors[:, order]
    Φtilde = C * U

    centre_table = reshape(centres_all, n.N_b, Np)
    r_translation = 0.0
    Lp = _scaled_physics(p, s).Lp
    for q = 1:(Np-1), a = 1:n.N_b
        r_translation =
            max(r_translation, abs(centre_table[a, q+1] - centre_table[a, q] - Lp) / Lp)
    end
    qcentral = n.P_basis + 1
    r_translation_nearest = 0.0
    if qcentral < Np
        for a = 1:n.N_b
            r_translation_nearest = max(
                r_translation_nearest,
                abs(centre_table[a, qcentral+1] - centre_table[a, qcentral] - Lp) / Lp,
            )
        end
    end
    r_translation_two_pairs = r_translation_nearest
    if qcentral + 1 < Np
        for a = 1:n.N_b
            r_translation_two_pairs = max(
                r_translation_two_pairs,
                abs(centre_table[a, qcentral+2] - centre_table[a, qcentral+1] - Lp) / Lp,
            )
        end
    end

    rows = (n.P_basis*n.N_z+1):((n.P_basis+1)*n.N_z)
    columns = (n.P_basis*n.N_b+1):((n.P_basis+1)*n.N_b)
    B = Matrix(Φtilde[rows, columns])
    S = Hermitian((B' * B + (B' * B)') / 2)
    eigS = eigen(S)
    smin, smax = extrema(eigS.values)
    κS = smax / smin
    if smin ≤ 0 || smin ≤ 100eps(Float64) * smax || κS > 1e8
        throw(ArgumentError("cell overlap is singular or ill-conditioned: κ₂=$κS"))
    end
    Sinvhalf = eigS.vectors * Diagonal(1 ./ sqrt.(eigS.values)) * eigS.vectors'
    Φ = ComplexF64.(B * Sinvhalf)
    _phase_fix!(Φ)
    r_orth = norm(Φ' * Φ - I)

    # Use the central diagonal block of the *same* window operator.  Its two
    # edge diagonals therefore contain the harmonic seam hopping.  Rebuilding
    # an isolated one-period Dirichlet matrix would use the wrong face mass.
    Hcell = Matrix(Hwin[rows, rows])
    H₀ = Matrix(Φ' * Hcell * Φ)
    Z = Matrix(Φ' * Diagonal(grids.x) * Φ)
    M⁻¹ = Matrix(Φ' * Diagonal(1 ./ profiles.m_parallelᵣ) * Φ)
    right_rows = ((n.P_basis+1)*n.N_z+1):((n.P_basis+2)*n.N_z)
    Hright = Matrix(Hwin[rows, right_rows])
    T₊ = Matrix(Φ' * Hright * Φ)
    T₋ = Matrix(T₊')

    χ = Φ ./ sqrt(Δx)
    centres = real.(diag(Z))
    second = real.(diag(Φ' * Diagonal(grids.x .^ 2) * Φ))
    spreads = sqrt.(max.(second .- centres .^ 2, 0.0))
    return BasisData(
        localization,
        Φ,
        χ,
        centres,
        spreads,
        collect(eigS.values),
        κS,
        H₀,
        Z,
        M⁻¹,
        Matrix(T₊),
        T₋,
        collect(eigenvalues_window),
        Egrid₀,
        r_orth,
        r_eigen,
        r_translation,
        r_translation_nearest,
        r_translation_two_pairs,
    )
end

"""Build the configured basis through a stable numerical-strategy interface.

`localization=:pzp` is the truncated
projected-position/Löwdin approximation and is not a low-rank physical oracle;
`localization=:none` uses central-cell energy eigenstates and is a controlled
numerical choice whenever the basis is truncated. `:real_space` preserves the
entire cell FD space and requires `N_b == N_z`.
"""
build_basis(
    p::PhysicalParameters,
    n::NumericalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    profiles::MaterialProfiles;
    localization::Symbol = :pzp,
) = build_localized_basis(p, n, s, grids, profiles; localization)

"""
    project_hamiltonians(problem, Uᴴ)

Project the periodic Hartree energy `Uᴴ(x_g)` and assemble the full dense
matrix `h(κ_m;Uᴴ)` for every radial momentum.  The returned array order is
`(N_k,N_b,N_b)` and every value is dimensionless.

Implements [EQ-HRED-001](@ref eq-reduced-hamiltonian).
See [BenDaniel–Duke and the fixed PzP basis](@ref theory-basis) for every
term of the projected Hamiltonian.
"""
function project_hamiltonians(
    basis::BasisData,
    grids::ModelGrids,
    profiles::MaterialProfiles,
    p::PhysicalParameters,
    s::ScaleSystem,
    Uᴴ::AbstractVector{<:Real},
)
    length(Uᴴ) == length(grids.x) || throw(DimensionMismatch("Uᴴ must have N_z entries"))
    Nb = size(basis.Φ, 2)
    Nk = length(grids.κ)
    sp = _scaled_physics(p, s)
    Uproj = basis.Φ' * Diagonal(Uᴴ) * basis.Φ
    Zshifted = basis.Z - sp.z₀ * Matrix{ComplexF64}(I, Nb, Nb)
    h = zeros(ComplexF64, Nk, Nb, Nb)
    for m = 1:Nk
        hm = basis.H₀ - sp.F * Zshifted + Uproj + _C̄(s) * grids.κ[m]^2 * basis.M⁻¹
        h[m, :, :] .= hm
    end
    return h
end

function project_hamiltonians(problem::NEGFProblem, Uᴴ::AbstractVector{<:Real})
    if problem.models.dispersion === :parabolic
        return project_hamiltonians(
            problem.basis,
            problem.grids,
            problem.profiles,
            problem.physical,
            problem.scales,
            Uᴴ,
        )
    end
    length(Uᴴ) == problem.numerical.N_z && all(isfinite, Uᴴ) ||
        throw(ArgumentError("Hartree field must be finite and match the spatial grid"))
    basis = problem.basis
    nb = size(basis.Φ, 2)
    sp = _scaled_physics(problem.physical, problem.scales)
    static =
        basis.H₀-sp.F*(basis.Z-sp.z₀*Matrix{ComplexF64}(I, nb, nb)) +
        basis.Φ'*Diagonal(Uᴴ)*basis.Φ
    h = zeros(ComplexF64, length(problem.grids.κ), nb, nb)
    for m in eachindex(problem.grids.κ)
        h[m, :, :] .= static+_transverse_hamiltonian(problem, problem.grids.κ[m])
    end
    return h
end

"""
    build_multiperiod_reference_basis(p, n, scales, grids, profiles;
                                      localization=:pzp_tails)

Small, complete-window reference for comparing `:pzp_tails`,
`:bloch_wannier`, and `:wannier_stark`. The window has `2P_basis+1` periods.
All envelope tails and the full projected Hamiltonian are retained, without
forcing a nearest-neighbour or rank-one coupling. `transform` acts on the
entire real-space window. This return type is intentionally distinct from
cell-local `BasisData`: it cannot accidentally be used with cell-local
density/scattering/embedding operators. Use `project_reference_operator`
to transform *every* full-window observable/vertex consistently.

PzP uses a zero-field Dirichlet spectral subspace. Wannier uses the low bands
of the periodic BDD cell, parallel-transported with a distributed unitary
holonomy before a discrete Bloch transform. Wannier--Stark uses the actual
finite-window field Hamiltonian; finite-window edge states and truncation
still require a domain-size and energy-window convergence study.
"""
function build_multiperiod_reference_basis(
    p::PhysicalParameters,
    n::NumericalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    profiles::MaterialProfiles;
    localization::Symbol = :pzp_tails,
)
    localization in (:pzp_tails, :bloch_wannier, :wannier_stark) ||
        throw(ArgumentError("unsupported full-window localization $localization"))
    count, Nz, Nb = 2n.P_basis + 1, n.N_z, n.N_b
    count >= 3 && 1 <= Nb <= Nz ||
        throw(ArgumentError("invalid reference window dimensions"))
    dx = only(unique(grids.wˣ))
    sp = _scaled_physics(p, s)
    coordinates = [(r - n.P_basis) * sp.Lp + grids.x[g] for r = 0:(count-1) for g = 1:Nz]
    boundary = localization === :bloch_wannier ? :periodic : :dirichlet
    Hfd = build_bdd_hamiltonian(
        repeat(profiles.Eᶜ, count),
        repeat(profiles.mᶻᵣ, count),
        dx,
        _C̄(s);
        boundary,
    )
    if localization === :wannier_stark
        Hfd .-= Diagonal(sp.F .* (coordinates .- sp.z₀))
    end
    Nstates = count * Nb
    gauge_minimum_overlap = nothing
    if localization === :bloch_wannier
        gauge_minimum_overlap = 1.0
        rows = (n.P_basis*Nz+1):((n.P_basis+1)*Nz)
        next = ((n.P_basis+1)*Nz+1):((n.P_basis+2)*Nz)
        Hcell, T = Hfd[rows, rows], Hfd[rows, next]
        phases = [2π*j/count for j = 0:(count-1)]
        frames = Matrix{ComplexF64}[]
        band_energies = Vector{Float64}[]
        for q in phases
            eig = eigen(Hermitian(Hcell + T * exp(im*q) + T' * exp(-im*q)))
            frame = ComplexF64.(eig.vectors[:, 1:Nb])
            push!(band_energies, eig.values[1:Nb])
            if !isempty(frames)
                polar = svd(frames[end]' * frame)
                gauge_minimum_overlap = min(gauge_minimum_overlap, minimum(polar.S))
                frame = frame * polar.V * polar.U'
            end
            push!(frames, frame)
        end
        closure = svd(frames[end]' * frames[1])
        gauge_minimum_overlap = min(gauge_minimum_overlap, minimum(closure.S))
        holonomy = closure.V * closure.U'
        phase_eigen = eigen(holonomy)
        for j = 1:count
            correction =
                phase_eigen.vectors *
                Diagonal(exp.(-im .* angle.(phase_eigen.values) .* ((j-1)/count))) *
                inv(phase_eigen.vectors)
            frames[j] = frames[j] * correction
        end
        Φ = zeros(ComplexF64, count*Nz, Nstates)
        for r = 0:(count-1), R = 0:(count-1), j = 1:count
            rr, cc = (r*Nz+1):((r+1)*Nz), (R*Nb+1):((R+1)*Nb)
            Φ[rr, cc] .+= exp(im * phases[j] * (r-R)) .* frames[j] ./ count
        end
        reference_energies = sort(vcat(band_energies...))
    else
        eig = eigen(Hermitian(Hfd))
        Φ = ComplexF64.(eig.vectors[:, 1:Nstates])
        reference_energies = eig.values[1:Nstates]
        if localization === :pzp_tails
            projected_position = Φ' * Diagonal(coordinates) * Φ
            localized = eigen(Hermitian(projected_position))
            Φ = Φ * localized.vectors
        else
            # Energy eigenfunctions are field-adapted; order only for a
            # meaningful central-window selection, not a second localisation.
            centres = real.(diag(Φ' * Diagonal(coordinates) * Φ))
            Φ = Φ[:, sortperm(centres)]
        end
    end
    _phase_fix!(Φ)
    H = Φ' * Hfd * Φ
    Z = Φ' * Diagonal(coordinates) * Φ
    centres = real.(diag(Z))
    spreads =
        sqrt.(max.(real.(diag(Φ' * Diagonal(coordinates .^ 2) * Φ)) - centres .^ 2, 0.0))
    central_rows = (n.P_basis*Nz+1):((n.P_basis+1)*Nz)
    central_columns = (n.P_basis*Nb+1):((n.P_basis+1)*Nb)
    tail_weights = 1 .- vec(sum(abs2, Φ[central_rows, central_columns]; dims = 1))
    residual = norm(Hfd * Φ - Φ * H) / max(norm(Hfd * Φ), eps(Float64))
    return (
        localization = localization,
        periods = count,
        boundary = boundary,
        transform = Φ,
        envelopes = Φ ./ sqrt(dx),
        coordinates = coordinates,
        weights = fill(dx, length(coordinates)),
        hamiltonian = H,
        position = Z,
        inverse_parallel_mass = Φ' * Diagonal(repeat(1 ./ profiles.m_parallelᵣ, count)) * Φ,
        real_space_hamiltonian = Hfd,
        centres = centres,
        spreads = spreads,
        central_rows = central_rows,
        central_columns = central_columns,
        central_tail_weights = tail_weights,
        overlap = Φ' * Φ,
        subspace_residual = residual,
        reference_energies = reference_energies,
        field_included = localization === :wannier_stark,
        operator_scope = :full_window,
        coupling_truncation = :none,
        gauge_minimum_overlap = gauge_minimum_overlap,
        gauge_transport_regular = gauge_minimum_overlap === nothing ? nothing :
                                  gauge_minimum_overlap > sqrt(eps(Float64)),
    )
end

"""Project a full-window Hamiltonian, density/current operator or scattering vertex."""
function project_reference_operator(reference, operator::AbstractMatrix)
    N = size(reference.transform, 1)
    size(operator) == (N, N) ||
        throw(DimensionMismatch("operator must cover the entire reference window"))
    return reference.transform' * operator * reference.transform
end

"""Measure low-band errors of a cell-local reduction against its original FD Bloch operator."""
function basis_dispersion_diagnostic(
    basis::BasisData,
    profiles::MaterialProfiles,
    grids::ModelGrids,
    scales::ScaleSystem;
    phases = range(-π, π; length = 17),
)
    full = build_bdd_hamiltonian(profiles, grids, scales; boundary = :periodic)
    seam = _boundary_hopping(profiles, grids, scales)
    full[1, end] = full[end, 1] = 0
    T = zeros(ComplexF64, size(full))
    T[end, 1] = -seam
    Nb = size(basis.Φ, 2)
    errors = zeros(Float64, length(phases), Nb)
    for (j, q) in enumerate(phases)
        exact = eigvals(Hermitian(full + exp(im*q) * T + exp(-im*q) * T'))[1:Nb]
        reduced =
            eigvals(Hermitian(basis.H₀ + exp(im*q) * basis.T₊ + exp(-im*q) * basis.T₋))
        errors[j, :] .= reduced - exact
    end
    return (
        phases = collect(phases),
        energy_error = errors,
        maximum_error_eV = maximum(abs, errors) * scales.E₀_eV,
        comparison = :same_discrete_bdd_operator,
    )
end

"""
Cheap pre-SCBA spectrum/window check. Local k-resolved poles are necessary
diagnostics, not the full interacting spectrum. The coupling norm supplies
an additional zero-field nearest-neighbour spectral bound. No small
population argument removes high unoccupied poles from a matrix sum rule.
"""
function spectral_window_diagnostic(problem::NEGFProblem, h::AbstractArray{<:Number,3})
    size(h, 1) == length(problem.grids.κ) ||
        throw(DimensionMismatch("momentum axis differs"))
    spectra = [eigvals(Hermitian(Matrix(view(h, m, :, :)))) for m in axes(h, 1)]
    minimum_level = minimum(minimum, spectra)
    maximum_level = maximum(maximum, spectra)
    coupling_margin = opnorm(problem.basis.T₊) + opnorm(problem.basis.T₋)
    emin, emax = first(problem.grids.ε), last(problem.grids.ε)
    trusted = problem.grids.ε[problem.grids.trusted_energy]
    return (
        local_minimum = minimum_level,
        local_maximum = maximum_level,
        sampled_minimum = emin,
        sampled_maximum = emax,
        trusted_minimum = first(trusted),
        trusted_maximum = last(trusted),
        local_spectrum_covered = emin <= minimum_level && maximum_level <= emax,
        trusted_local_spectrum_covered = first(trusted) <= minimum_level &&
                                         maximum_level <= last(trusted),
        coupling_margin = coupling_margin,
        estimate = :local_poles_not_interacting_spectrum,
    )
end
