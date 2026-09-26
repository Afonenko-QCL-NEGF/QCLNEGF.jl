"""
    OpticalResponse

Frequency-domain, small-signal optical response of one converged stationary
NEGF state.  `photon_energy` is in eV, `frequency` in Hz,
`angular_frequency` in s^-1, `gain` in m^-1 and `susceptibility` and
`refractive_index` are dimensionless.  Positive `gain` denotes amplification;
negative `gain` denotes absorption.

The current implementation is the length-gauge *bare bubble*.  It retains
the complete matrix Green functions and their SCBA broadening, but deliberately
sets the optical variation of every self-energy to zero.  Consequently
`vertex_corrections == false`; this fact is part of the returned data rather
than an undocumented approximation.

See [Optical response](@ref theory-optical-response),
[Production backend](@ref theory-production), and
[debug/final visualization](@ref theory-visualization).
"""
struct OpticalResponse
    photon_energy::Vector{EnergyQuantity}
    frequency::Vector{typeof(1.0u"Hz")}
    angular_frequency::Vector{typeof(1.0u"s^-1")}
    susceptibility::Vector{ComplexF64}
    refractive_index::Vector{ComplexF64}
    gain::Vector{typeof(1.0u"m^-1")}
    trusted::BitVector
    edge_loss::Vector{Float64}
    background_permittivity::Float64
    model::Symbol
    gauge::Symbol
    vertex_corrections::Bool
end

function Base.show(io::IO, response::OpticalResponse)
    trusted = count(response.trusted)
    print(
        io,
        "OpticalResponse(model=$(response.model), ",
        "N_photon=$(length(response.photon_energy)), trusted=$trusted)",
    )
end

"""
    _optical_shift_plan(ε, δ)

Return `(left, fraction, valid)` for the non-periodic piecewise-linear
evaluation `X(ε[e] + δ)`.  Every valid entry is represented as
`(1-fraction[e])X[left[e]] + fraction[e]X[left[e]+1]`; values outside the
energy interval are invalid and are never wrapped.  This O(N_E) plan is the
optical analogue of [`build_shift_matrix`](@ref), without allocating a dense
`N_E × N_E` matrix.
"""
function _optical_shift_plan(ε::AbstractVector{<:Real}, δ::Real)
    issorted(ε) || throw(ArgumentError("energy nodes must be sorted"))
    δ ≥ 0 || throw(ArgumentError("photon-energy shift must be nonnegative"))
    N = length(ε)
    N ≥ 2 || throw(ArgumentError("at least two energy nodes are required"))
    left = zeros(Int, N)
    fraction = zeros(Float64, N)
    valid = falses(N)
    for e in eachindex(ε)
        target = ε[e] + δ
        (target < ε[1] || target > ε[end]) && continue
        j = target == ε[end] ? N - 1 : clamp(searchsortedlast(ε, target), 1, N - 1)
        θ = (target - ε[j]) / (ε[j+1] - ε[j])
        left[e] = j
        fraction[e] = θ
        valid[e] = true
    end
    return left, fraction, valid
end

@inline function _trace_product(A::AbstractMatrix, B::AbstractMatrix)
    size(A) == reverse(size(B)) || throw(DimensionMismatch("trace-product matrices differ"))
    return _trace_product_unchecked(A, B)
end

@inline function _trace_product_unchecked(A::AbstractMatrix, B::AbstractMatrix)
    value = zero(promote_type(eltype(A), eltype(B)))
    @inbounds for a in axes(A, 1), b in axes(A, 2)
        value += A[a, b] * B[b, a]
    end
    return value
end

function _interpolate_block!(
    destination::AbstractMatrix,
    X::AbstractArray{<:Number,4},
    j::Int,
    m::Int,
    θ::Real,
)
    Nb = size(destination, 1)
    size(destination, 2) == Nb || throw(DimensionMismatch("destination must be square"))
    return _interpolate_block_unchecked!(destination, X, j, m, θ)
end

@inline function _interpolate_block_unchecked!(
    destination::AbstractMatrix,
    X::AbstractArray{<:Number,4},
    j::Int,
    m::Int,
    θ::Real,
)
    Nb = size(destination, 1)
    @inbounds for a = 1:Nb, b = 1:Nb
        destination[a, b] = (1 - θ) * X[j, m, a, b] + θ * X[j+1, m, a, b]
    end
    return destination
end

@doc raw"""
    _bare_bubble_sum(ε, wᴱ, wᵏ, Z, Gᴿ, Gˡ, δ)

Evaluate the dimensionless Green-function part of the length-gauge bare
bubble at positive dimensionless photon energy `δ`.  The returned tuple is
`(S, valid)`, where

```math
S(\delta)=\sum_{e,m}w^E_e w^k_m\operatorname{tr}\!\left\{
Z\left[G^R(E_e+\delta)ZG^<(E_e)
+G^<(E_e+\delta)ZG^A(E_e)\right]\right\}.
```

`valid[e]` records whether `E_e+δ` lies inside the stored energy interval.
The routine performs no cyclic extrapolation and allocates only four
`N_b × N_b` work matrices.
"""
function _bare_bubble_sum(
    ε::AbstractVector{<:Real},
    wᴱ::AbstractVector{<:Real},
    wᵏ::AbstractVector{<:Real},
    Z::AbstractMatrix{<:Number},
    Gᴿ::AbstractArray{<:Number,4},
    Gˡ::AbstractArray{<:Number,4},
    δ::Real,
)
    size(Gᴿ) == size(Gˡ) || throw(DimensionMismatch("Gᴿ and Gˡ differ"))
    NE, Nk, Nb, Nb2 = size(Gᴿ)
    Nb == Nb2 || throw(DimensionMismatch("Green blocks must be square"))
    length(ε) == NE == length(wᴱ) || throw(DimensionMismatch("energy axes differ"))
    length(wᵏ) == Nk || throw(DimensionMismatch("momentum axes differ"))
    size(Z) == (Nb, Nb) || throw(DimensionMismatch("dipole matrix has wrong shape"))
    left, fraction, valid = _optical_shift_plan(ε, δ)

    Gᴿplus = Matrix{ComplexF64}(undef, Nb, Nb)
    Gˡplus = Matrix{ComplexF64}(undef, Nb, Nb)
    work₁ = Matrix{ComplexF64}(undef, Nb, Nb)
    work₂ = Matrix{ComplexF64}(undef, Nb, Nb)
    Zc = ComplexF64.(Z)
    response = 0.0 + 0.0im
    for e = 1:NE
        valid[e] || continue
        j = left[e]
        θ = fraction[e]
        for m = 1:Nk
            _interpolate_block!(Gᴿplus, Gᴿ, j, m, θ)
            _interpolate_block!(Gˡplus, Gˡ, j, m, θ)

            # Gᴿ(E+ℏω) Z G<(E)
            mul!(work₁, Gᴿplus, Zc)
            mul!(work₂, work₁, view(Gˡ, e, m, :, :))
            term = _trace_product(Zc, work₂)

            # G<(E+ℏω) Z Gᴬ(E), with Gᴬ=(Gᴿ)†.
            mul!(work₁, Gˡplus, Zc)
            mul!(work₂, work₁, adjoint(view(Gᴿ, e, m, :, :)))
            term += _trace_product(Zc, work₂)
            response += wᴱ[e] * wᵏ[m] * term
        end
    end
    return response, valid
end

function _dipole_sandwich_blocks(Z::AbstractMatrix{<:Number}, X::AbstractArray{<:Number,4})
    NE, Nk, Nb, Nb2 = size(X)
    Nb == Nb2 || throw(DimensionMismatch("matrix blocks must be square"))
    size(Z) == (Nb, Nb) || throw(DimensionMismatch("dipole matrix has wrong shape"))
    Zc = ComplexF64.(Z)
    Y = zeros(ComplexF64, size(X))
    work₁ = Matrix{ComplexF64}(undef, Nb, Nb)
    work₂ = similar(work₁)
    for e = 1:NE, m = 1:Nk
        mul!(work₁, Zc, view(X, e, m, :, :))
        mul!(work₂, work₁, Zc)
        view(Y, e, m, :, :) .= work₂
    end
    return Y
end

function _bare_bubble_sum_cached(
    ε::AbstractVector{<:Real},
    wᴱ::AbstractVector{<:Real},
    wᵏ::AbstractVector{<:Real},
    ZGᴿZ::AbstractArray{<:Number,4},
    ZGˡZ::AbstractArray{<:Number,4},
    Gᴿ::AbstractArray{<:Number,4},
    Gˡ::AbstractArray{<:Number,4},
    δ::Real,
)
    size(Gᴿ) == size(Gˡ) == size(ZGᴿZ) == size(ZGˡZ) ||
        throw(DimensionMismatch("cached optical arrays differ"))
    NE, Nk, Nb, Nb2 = size(Gᴿ)
    Nb == Nb2 || throw(DimensionMismatch("Green blocks must be square"))
    length(ε) == NE == length(wᴱ) || throw(DimensionMismatch("energy axes differ"))
    length(wᵏ) == Nk || throw(DimensionMismatch("momentum axes differ"))
    left, fraction, valid = _optical_shift_plan(ε, δ)
    ZGᴿZplus = Matrix{ComplexF64}(undef, Nb, Nb)
    ZGˡZplus = similar(ZGᴿZplus)
    response = 0.0 + 0.0im
    for e = 1:NE
        valid[e] || continue
        j = left[e]
        θ = fraction[e]
        for m = 1:Nk
            _interpolate_block_unchecked!(ZGᴿZplus, ZGᴿZ, j, m, θ)
            _interpolate_block_unchecked!(ZGˡZplus, ZGˡZ, j, m, θ)
            term = _trace_product_unchecked(ZGᴿZplus, view(Gˡ, e, m, :, :))
            term += _trace_product_unchecked(ZGˡZplus, adjoint(view(Gᴿ, e, m, :, :)))
            response += wᴱ[e] * wᵏ[m] * term
        end
    end
    return response, valid
end

function _spectral_weight_by_energy(problem::NEGFProblem, green::GreenState)
    NE, Nk, Nb, _ = size(green.A)
    weight = zeros(Float64, NE)
    for e = 1:NE, m = 1:Nk
        # `abs` makes the edge diagnostic well defined even before the final
        # PSD tolerance is met.  It is a diagnostic only and never alters χ.
        weight[e] += problem.grids.wᵏ[m] * abs(real(tr(view(green.A, e, m, :, :))))
    end
    return weight
end

@doc raw"""
    bare_bubble_optical_response(problem, green, photon_energies;
        background_permittivity=problem.physical.ε_s,
        edge_tolerance=1e-4, threaded=true)

Calculate the z-polarized small-signal material gain of a converged
stationary state.  `photon_energies` must contain positive Unitful energies.
The result uses the full matrix `Gᴿ(E,k)` and `Gˡ(E,k)` and the projected
intracell position matrix `Z`; no subband-population or Lorentzian-linewidth
ansatz is introduced.

For a reference field `F_ref = 1 V/m`, the dimensionless perturbation is

```math
\delta \widehat U
= \frac{eF_{ref}L_0}{E_0}\,Z,
```

and the bare-bubble change is

```math
\delta G^<(E;\omega)=G^R(E+\hbar\omega)\delta\widehat U G^<(E)
+G^<(E+\hbar\omega)\delta\widehat U G^A(E).
```

The sheet density response is reconstructed with exactly the same
`wᴱ`, `wᵏ` and spin convention as [`sheet_density_matrix`](@ref).  With
electron charge `-e`, one-period polarization and susceptibility are

```math
\delta P_z=-\frac{e}{L_pL_0}\operatorname{tr}(Z\delta\bar\rho),
\qquad \chi=\frac{\delta P_z}{\epsilon_0F_{ref}}.
```

Finally, for the `exp(-iωt)` convention,

```math
\widetilde n=\sqrt{\epsilon_b+\chi},\qquad
g(\omega)=-\frac{2\omega}{c}\operatorname{Im}\widetilde n.
```

Thus positive `g` means gain and negative `g` means absorption.  The exact
square-root expression is used rather than its weak-gain expansion.

This is **not** the gauge-invariant conserving response of the complete SCBA
model: optical self-energy variations (ladder/vertex corrections),
neighbor-period dipole blocks, photon-assisted changes of Poisson, and cavity
mode confinement are omitted.  It is suitable for a reproducible diagnostic
and for gain-frequency trends, but must not be presented as a quantitative
reproduction of the QCL gain without convergence and experimental
calibration.

`trusted[q]` is true only if photon energy `q` is no larger than the configured
Hilbert trust margin and a symmetric source/target spectral-weight proxy for
out-of-window or Hilbert-edge-contaminated `(E,E+ℏω)` pairs is at most
`edge_tolerance`.  This proxy samples both ``A(E)`` and ``A(E+\hbar\omega)``;
it is deliberately conservative but is not a rigorous bound on the omitted
optical integrand.  All points are returned so edge sensitivity cannot be
silently hidden.

See [Optical response](@ref theory-optical-response) for the response
equations and limitations, and [Production backend](@ref theory-production)
for its role in the reference design workflow.
"""
function bare_bubble_optical_response(
    problem::NEGFProblem,
    green::GreenState,
    photon_energies;
    background_permittivity::Real = problem.physical.ε_s,
    edge_tolerance::Real = 1e-4,
    threaded::Bool = true,
)
    background_permittivity > 0 ||
        throw(ArgumentError("background_permittivity must be positive"))
    0 ≤ edge_tolerance < 1 || throw(ArgumentError("edge_tolerance must lie in [0,1)"))
    expected = (
        problem.numerical.N_E,
        problem.numerical.N_k,
        problem.numerical.N_b,
        problem.numerical.N_b,
    )
    size(green.Gᴿ) == expected == size(green.Gˡ) == size(green.A) ||
        throw(DimensionMismatch("Green state does not match the problem"))

    energies = EnergyQuantity[_energy(value) for value in photon_energies]
    isempty(energies) && throw(ArgumentError("photon-energy grid is empty"))
    all(value -> value > 0u"eV", energies) ||
        throw(ArgumentError("all photon energies must be positive"))
    energy_values = _electronvolts.(energies)
    all(isfinite, energy_values) || throw(ArgumentError("photon energies must be finite"))
    all(diff(energy_values) .> 0) ||
        throw(ArgumentError("photon energies must be strictly increasing"))

    Nb = problem.numerical.N_b
    Z = ComplexF64.(problem.basis.Z)
    rZ = norm(Z - Z') / (norm(Z) + eps(Float64))
    rZ ≤ _OPTICAL_POSITION_HERMITICITY_LIMIT ||
        throw(DomainError(rZ, "projected position matrix is not Hermitian"))
    # Fix the arbitrary intracell coordinate origin.  In a complete conserving
    # response this identity shift cancels analytically; centering also limits
    # finite-window roundoff in the bare bubble.
    Z .-= real(tr(Z)) / Nb .* Matrix{ComplexF64}(I, Nb, Nb)

    Fref = 1.0u"V/m"
    drive_scale = Float64(
        ustrip(
            Unitful.NoUnits,
            uconvert(
                Unitful.NoUnits,
                CODATA.e * Fref * problem.scales.L₀ / uconvert(u"J", problem.scales.E₀),
            ),
        ),
    )
    Lp = period_length(problem.physical)
    polarization_scale = Float64(
        ustrip(
            Unitful.NoUnits,
            uconvert(
                Unitful.NoUnits,
                -CODATA.e / (Lp * problem.scales.L₀) / (CODATA.ε₀ * Fref),
            ),
        ),
    )

    Nω = length(energies)
    χ = zeros(ComplexF64, Nω)
    edge_loss = zeros(Float64, Nω)
    trusted = falses(Nω)
    spectral_by_energy = _spectral_weight_by_energy(problem, green)
    spectral_total = dot(problem.grids.wᴱ, spectral_by_energy)
    spectral_total > 0 ||
        throw(DomainError(spectral_total, "spectral weight must be positive"))
    margin_eV = _electronvolts(problem.numerical.M_E)

    # Two cached matrix fields change the photon sweep from O(Nω NE Nk Nb³)
    # to O(NE Nk Nb³ + Nω NE Nk Nb²).  At the reference design baseline each cache is
    # about 56 MiB; this deterministic 112 MiB trade-off is explicitly bounded
    # and avoids allocating shifted four-dimensional fields for every photon.
    ZGᴿZ = _dipole_sandwich_blocks(Z, green.Gᴿ)
    ZGˡZ = _dipole_sandwich_blocks(Z, green.Gˡ)

    function evaluate!(q)
        δeV = _electronvolts(energies[q])
        δ = δeV / problem.scales.E₀_eV
        bubble, valid = _bare_bubble_sum_cached(
            problem.grids.ε,
            problem.grids.wᴱ,
            problem.grids.wᵏ,
            ZGᴿZ,
            ZGˡZ,
            green.Gᴿ,
            green.Gˡ,
            δ,
        )
        density_response = (-im * problem.physical.g_s / (2π)) * drive_scale * bubble
        χ[q] = polarization_scale * density_response
        left, fraction, _ = _optical_shift_plan(problem.grids.ε, δ)
        missing_proxy = 0.0
        total_proxy = 0.0
        for e in eachindex(valid)
            pair_trusted = valid[e] && problem.grids.trusted_energy[e]
            target_weight = spectral_by_energy[end]
            if pair_trusted
                j = left[e]
                θ = fraction[e]
                pair_trusted &=
                    (θ == 1 || problem.grids.trusted_energy[j]) &&
                    (θ == 0 || problem.grids.trusted_energy[j+1])
            end
            if valid[e]
                j = left[e]
                θ = fraction[e]
                target_weight =
                    (1 - θ) * spectral_by_energy[j] + θ * spectral_by_energy[j+1]
            end
            pair_proxy = problem.grids.wᴱ[e] * (spectral_by_energy[e] + target_weight)
            total_proxy += pair_proxy
            pair_trusted || (missing_proxy += pair_proxy)
        end
        edge_loss[q] = missing_proxy / max(total_proxy, eps(Float64))
        trusted[q] = δeV ≤ margin_eV + 10eps(margin_eV) && edge_loss[q] ≤ edge_tolerance
        return nothing
    end

    if threaded && Threads.nthreads() > 1 && Nω > 1
        Threads.@threads :static for q = 1:Nω
            evaluate!(q)
        end
    else
        for q = 1:Nω
            evaluate!(q)
        end
    end
    all(
        q -> isfinite(real(χ[q])) && isfinite(imag(χ[q])) && isfinite(edge_loss[q]),
        eachindex(χ),
    ) || throw(
        DomainError(
            (χ = χ, edge_loss = edge_loss),
            "optical response contains a non-finite value",
        ),
    )

    ntilde = sqrt.(complex(background_permittivity) .+ χ)
    angular_frequency = typeof(1.0u"s^-1")[]
    frequency = typeof(1.0u"Hz")[]
    gain = typeof(1.0u"m^-1")[]
    for q = 1:Nω
        ω = uconvert(u"s^-1", energies[q] / CODATA.ħₑᵥ)
        ν = uconvert(u"Hz", ω / (2π))
        g = uconvert(u"m^-1", -2imag(ntilde[q]) * ω / CODATA.c)
        push!(angular_frequency, ω)
        push!(frequency, ν)
        push!(gain, g)
    end
    return OpticalResponse(
        energies,
        frequency,
        angular_frequency,
        χ,
        ntilde,
        gain,
        trusted,
        edge_loss,
        Float64(background_permittivity),
        :bare_bubble,
        :length,
        false,
    )
end

"""
    bare_bubble_optical_response(solution, photon_energies;
                                 require_converged=true, kwargs...)

Convenience overload for a closed Poisson--SCBA solution.  Production plots
reject an unconverged or validation-failed state by default; set
`require_converged=false` only for explicitly labelled diagnostics.

See [Optical response](@ref theory-optical-response) and
[Production backend](@ref theory-production).
"""
function bare_bubble_optical_response(
    solution::NEGFSolution,
    photon_energies;
    require_converged::Bool = true,
    kwargs...,
)
    if require_converged && !solution.converged
        throw(
            ArgumentError(
                "optical response requires a converged, validated solution; " *
                "received status=$(solution.status)",
            ),
        )
    end
    return bare_bubble_optical_response(
        solution.problem,
        solution.scba.green,
        photon_energies;
        kwargs...,
    )
end

"""
    peak_gain(response; trusted_only=true)

Return the largest material gain and its photon energy/frequency.  By default
only points certified by `response.trusted` are eligible.

See [Optical response](@ref theory-optical-response) and
[production gain visualization](@ref theory-production).
"""
function peak_gain(response::OpticalResponse; trusted_only::Bool = true)
    eligible = trusted_only ? findall(response.trusted) : collect(eachindex(response.gain))
    isempty(eligible) && throw(ArgumentError("the optical response has no eligible points"))
    values = Float64[ustrip(u"m^-1", response.gain[q]) for q in eligible]
    q = eligible[argmax(values)]
    return (
        gain = response.gain[q],
        photon_energy = response.photon_energy[q],
        frequency = response.frequency[q],
        index = q,
        trusted = response.trusted[q],
        edge_loss = response.edge_loss[q],
    )
end
