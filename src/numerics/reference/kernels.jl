_K₀(s::ScaleSystem) = s.E₀_J^2 * s.L₀_m^2
_electron_charge() = Float64(ustrip(u"C", CODATA.e))
_epsilon_zero() = Float64(ustrip(u"F/m", CODATA.ε₀))

function _angular_average_kernel(grids::ModelGrids, Nb::Int, kernel_at_q)
    Nk = length(grids.κ)
    K = zeros(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
    for m = 1:Nk, mp = 1:Nk
        accumulator = zeros(ComplexF64, Nb, Nb, Nb, Nb)
        for φ in grids.φ
            q = sqrt(
                max(grids.κ[m]^2 + grids.κ[mp]^2 - 2grids.κ[m] * grids.κ[mp] * cos(φ), 0.0),
            )
            accumulator .+= kernel_at_q(q)
        end
        K[m, mp, :, :, :, :] .= accumulator ./ length(grids.φ)
    end
    return K
end

"""
    adaptive_angular_average(evaluate, k, kprime; relative_tolerance=1e-5,
                             maximum_depth=14, initial_panels=8)

Average a rotationally symmetric tensor `K(q)` over the in-plane angle.
Symmetry reduces the integral to `[0,π]`; `φ=π*u²` resolves forward scattering
near zero without changing the Coulomb model. Recursive Simpson integration
uses positive quadrature weights, preserving covariance positivity. The
fine/coarse difference gives an *estimated*, not certified, error. A failed
depth limit is reported explicitly. Independent tolerance/angle refinement
is still required for a production model.
"""
function adaptive_angular_average(
    evaluate,
    k::Real,
    kprime::Real;
    relative_tolerance::Real = 1e-5,
    maximum_depth::Int = 14,
    initial_panels::Int = 8,
)
    k >= 0 && kprime >= 0 && isfinite(k) && isfinite(kprime) ||
        throw(ArgumentError("invalid radial momenta"))
    0 < relative_tolerance < 1 && maximum_depth >= 1 && initial_panels >= 1 ||
        throw(ArgumentError("invalid angular quadrature controls"))
    cache = Dict{Float64,Any}()
    function sample(u)
        return get!(cache, Float64(u)) do
            phi = π * u^2
            q = sqrt((k-kprime)^2 + 4k*kprime*sin(phi/2)^2)
            2u .* evaluate(q)
        end
    end
    failed = Ref(false)
    function integrate(a, b, fa, fm, fb, coarse, depth)
        middle = (a+b)/2
        fl, fr = sample((a+middle)/2), sample((middle+b)/2)
        left = (middle-a)/6 .* (fa .+ 4 .* fl .+ fm)
        right = (b-middle)/6 .* (fm .+ 4 .* fr .+ fb)
        fine = left .+ right
        error = norm(fine .- coarse)/15
        if error <= relative_tolerance * max(norm(fine), floatmin(Float64))
            return fine, error
        elseif depth == maximum_depth
            failed[] = true
            return fine, error
        end
        lv, le = integrate(a, middle, fa, fl, fm, left, depth+1)
        rv, re = integrate(middle, b, fm, fr, fb, right, depth+1)
        return lv .+ rv, le + re
    end
    total = similar(sample(0.0))
    fill!(total, 0)
    total_error = 0.0
    for panel = 0:(initial_panels-1)
        a, b = panel/initial_panels, (panel+1)/initial_panels
        fa, fm, fb = sample(a), sample((a+b)/2), sample(b)
        coarse = (b-a)/6 .* (fa .+ 4 .* fm .+ fb)
        value, error = integrate(a, b, fa, fm, fb, coarse, 1)
        total .+= value
        total_error += error
    end
    residual = total_error / max(norm(total), floatmin(Float64))
    return (
        value = total,
        estimated_relative_error = residual,
        evaluations = length(cache),
        accepted = !failed[] && residual <= relative_tolerance,
        construction = :adaptive_transformed_simpson,
    )
end

"""
    impurity_kernel(physical, scales, grids, profiles, basis)

Construct the complete screened ionized-donor covariance tensor.  Donor
coordinates and electron vertices use the same representative-period
cell-centred quadrature as Poisson.  The result is dimensionless and ordered
`(m,m′,a,c,d,b)`.

Implements [EQ-KIMP-001](@ref eq-impurity-kernel).
See [Microscopic kernels](@ref theory-kernels) for screening, donor averaging,
and the six-index ordering.
"""
function impurity_kernel(
    p::PhysicalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    profiles::MaterialProfiles,
    basis::BasisData;
    angular::Symbol = :legacy_uniform,
    angular_tolerance::Real = 1e-5,
)
    Nb = size(basis.χ, 2)
    eC = _electron_charge()
    ε0 = _epsilon_zero()
    qs = _inverse_metres(p.q_s)
    L0 = s.L₀_m
    x = grids.x
    weights_donor = grids.wˣ .* profiles.Nᴰ ./ L0^2
    function kernel_at_q(qbar)
        q = qbar / L0
        prefactor = eC^2 / (2ε0 * p.ε_s * (q + qs)) # J m²
        Ksi = zeros(ComplexF64, Nb, Nb, Nb, Nb)
        for gD in eachindex(x)
            V = zeros(ComplexF64, Nb, Nb)
            for a = 1:Nb, c = 1:Nb
                integral = 0.0 + 0.0im
                for g in eachindex(x)
                    integral +=
                        grids.wˣ[g] *
                        conj(basis.χ[g, a]) *
                        exp(-qbar * abs(x[g] - x[gD])) *
                        basis.χ[g, c]
                end
                V[a, c] = prefactor * integral
            end
            for a = 1:Nb, c = 1:Nb, d = 1:Nb, b = 1:Nb
                Ksi[a, c, d, b] += weights_donor[gD] * V[a, c] * conj(V[b, d])
            end
        end
        return Ksi ./ _K₀(s)
    end
    angular in (:legacy_uniform, :adaptive) ||
        throw(ArgumentError("unknown impurity angular quadrature"))
    angular === :legacy_uniform && return _angular_average_kernel(grids, Nb, kernel_at_q)
    Nk = length(grids.κ)
    K = zeros(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
    for m = 1:Nk, mp = 1:Nk
        result = adaptive_angular_average(
            kernel_at_q,
            grids.κ[m],
            grids.κ[mp];
            relative_tolerance = angular_tolerance,
        )
        result.accepted ||
            throw(ErrorException("impurity angular quadrature failed at ($m,$mp)"))
        K[m, mp, :, :, :, :] .= result.value
    end
    return K
end

function _interface_value(χ::AbstractVector, x::AbstractVector, xint::Real)
    if xint == 0
        return χ[1]
    end
    exact = findfirst(
        y -> isapprox(y, xint; rtol = 0, atol = 16eps(Float64) * max(abs(xint), 1)),
        x,
    )
    exact !== nothing && return χ[exact]
    gminus = searchsortedlast(x, xint)
    1 ≤ gminus < length(x) || throw(ArgumentError("interface lies outside cell centres"))
    θ = (xint - x[gminus]) / (x[gminus+1] - x[gminus])
    return (1 - θ) * χ[gminus] + θ * χ[gminus+1]
end

function _band_jump(p::PhysicalParameters, zint::LengthQuantity)
    z = _metres(zint)
    Lp = _metres(period_length(p))
    if iszero(z)
        return uconvert(u"J", p.layers[1].Eᶜ - p.layers[end].Eᶜ)
    end
    edges = cumsum(_metres(layer.d) for layer in p.layers)
    j = findfirst(edge -> isapprox(edge, z; rtol = 0, atol = 64eps(Float64) * Lp), edges)
    j === nothing && throw(ArgumentError("IFR coordinate is not a material interface"))
    j < length(p.layers) ||
        throw(ArgumentError("period end must be represented by seam z=0"))
    return uconvert(u"J", p.layers[j+1].Eᶜ - p.layers[j].Eᶜ)
end

"""
    interface_roughness_kernel(physical, scales, grids, basis)

Independent-interface Gaussian IFR tensor.  The periodic seam is counted once
from the right, and donor-segment edges are excluded by construction.

Implements [EQ-KIFR-001](@ref eq-ifr-kernel).
See [Microscopic kernels](@ref theory-kernels) for the Gaussian interface
spectrum and interface form factors.
"""
function interface_roughness_kernel(
    p::PhysicalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    basis::BasisData,
)
    Nb = size(basis.χ, 2)
    Δ = _metres(p.Δᴵᶠᴿ)
    Λ = _metres(p.Λᴵᶠᴿ)
    L0 = s.L₀_m
    D = Vector{Matrix{ComplexF64}}(undef, length(p.interfaces))
    for (j, zint) in pairs(p.interfaces)
        xint = _metres(zint) / L0
        χint = [_interface_value(view(basis.χ, :, a), grids.x, xint) for a = 1:Nb]
        ΔEc = Float64(ustrip(u"J", _band_jump(p, zint)))
        D[j] = [ΔEc * conj(χint[a]) * χint[c] / L0 for a = 1:Nb, c = 1:Nb]
    end
    function kernel_at_q(qbar)
        q = qbar / L0
        spectrum = π * Δ^2 * Λ^2 * exp(-q^2 * Λ^2 / 4)
        Ksi = zeros(ComplexF64, Nb, Nb, Nb, Nb)
        for Dj in D, a = 1:Nb, c = 1:Nb, d = 1:Nb, b = 1:Nb
            Ksi[a, c, d, b] += spectrum * Dj[a, c] * conj(Dj[b, d])
        end
        return Ksi ./ _K₀(s)
    end
    return _angular_average_kernel(grids, Nb, kernel_at_q)
end

"""
    acoustic_kernel(physical, scales, grids, basis)

Full four-index deformation-potential kernel in the elastic equipartition
limit, copied over `(m,m′)` because it is momentum independent.

Implements [EQ-KAC-001](@ref eq-acoustic-kernel).
See [Microscopic kernels](@ref theory-kernels) for the elastic
equipartition approximation.
"""
function acoustic_kernel(
    p::PhysicalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    basis::BasisData,
)
    Nb = size(basis.χ, 2)
    Nk = length(grids.κ)
    Ξ = Float64(ustrip(u"J", uconvert(u"J", p.Ξ)))
    kBT = Float64(ustrip(u"J", CODATA.kᴮ * p.Tᴸ))
    ρ = Float64(ustrip(u"kg/m^3", p.ρ_m))
    vs = Float64(ustrip(u"m/s", p.v_s))
    prefactor = Ξ^2 * kBT / (ρ * vs^2) # J² m³
    Kblock = zeros(ComplexF64, Nb, Nb, Nb, Nb)
    for a = 1:Nb, c = 1:Nb, d = 1:Nb, b = 1:Nb
        integral = 0.0 + 0.0im
        for g in eachindex(grids.x)
            integral +=
                grids.wˣ[g] *
                conj(basis.χ[g, a]) *
                basis.χ[g, c] *
                conj(basis.χ[g, d]) *
                basis.χ[g, b] / s.L₀_m
        end
        Kblock[a, c, d, b] = prefactor * integral / _K₀(s)
    end
    K = zeros(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
    for m = 1:Nk, mp = 1:Nk
        K[m, mp, :, :, :, :] .= Kblock
    end
    return K
end

"""
Full alloy-disorder tensor, or an exact zero tensor when the channel is
disabled. Implements [EQ-KALLOY-001](@ref eq-alloy-kernel).

See [Microscopic kernels](@ref theory-kernels) for the alloy covariance and
material-profile dependence.
"""
function alloy_kernel(
    p::PhysicalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    profiles::MaterialProfiles,
    basis::BasisData,
)
    Nb = size(basis.χ, 2)
    Nk = length(grids.κ)
    K = zeros(ComplexF64, Nk, Nk, Nb, Nb, Nb, Nb)
    (p.ΔV_alloy === nothing || p.Ω₀ === nothing) && return K
    ΔV = Float64(ustrip(u"J", uconvert(u"J", p.ΔV_alloy)))
    Ω = Float64(ustrip(u"m^3", p.Ω₀))
    block = zeros(ComplexF64, Nb, Nb, Nb, Nb)
    for a = 1:Nb, c = 1:Nb, d = 1:Nb, b = 1:Nb
        integral = 0.0 + 0.0im
        for g in eachindex(grids.x)
            alloy = profiles.x_Al[g] * (1 - profiles.x_Al[g])
            integral +=
                grids.wˣ[g] *
                alloy *
                conj(basis.χ[g, a]) *
                basis.χ[g, c] *
                conj(basis.χ[g, d]) *
                basis.χ[g, b] / s.L₀_m
        end
        block[a, c, d, b] = Ω * ΔV^2 * integral / _K₀(s)
    end
    for m = 1:Nk, mp = 1:Nk
        K[m, mp, :, :, :, :] .= block
    end
    return K
end

"""
    lo_phonon_kernel(physical, scales, grids, basis)

Screened bulk-like Fröhlich kernel and longitudinal form factors.  The direct
trapezoidal `q_z` rule from the theory is retained exactly.

Implements [EQ-KLO-001](@ref eq-lo-kernel).
See [Microscopic kernels](@ref theory-kernels) for the Fröhlich coupling,
phonon energy shifts, and longitudinal form factors.
"""
function lo_phonon_kernel(
    p::PhysicalParameters,
    s::ScaleSystem,
    grids::ModelGrids,
    basis::BasisData,
)
    Nb = size(basis.χ, 2)
    Nq = length(grids.qᶻ)
    Fᴸᴼ = zeros(ComplexF64, Nq, Nb, Nb)
    for h = 1:Nq, a = 1:Nb, c = 1:Nb
        for g in eachindex(grids.x)
            Fᴸᴼ[h, a, c] +=
                grids.wˣ[g] *
                conj(basis.χ[g, a]) *
                exp(im * grids.qᶻ[h] * grids.x[g]) *
                basis.χ[g, c]
        end
    end
    eC = _electron_charge()
    ε0 = _epsilon_zero()
    ħω = Float64(ustrip(u"J", uconvert(u"J", p.ħωᴸᴼ)))
    prefactor = eC^2 * ħω / (2ε0) * (1 / p.ε_∞ - 1 / p.ε_s) # J² m
    qscreen = scaled_value(p.qᴸᴼ_s, s, :wavenumber)
    function kernel_at_q(qbar)
        Ksi = zeros(ComplexF64, Nb, Nb, Nb, Nb)
        for h = 1:Nq
            denominator = qbar^2 + grids.qᶻ[h]^2 + qscreen^2
            factor = grids.wᑫᶻ[h] * s.L₀_m / (2π * denominator)
            for a = 1:Nb, c = 1:Nb, d = 1:Nb, b = 1:Nb
                Ksi[a, c, d, b] += factor * Fᴸᴼ[h, a, c] * conj(Fᴸᴼ[h, b, d])
            end
        end
        return prefactor .* Ksi ./ _K₀(s)
    end
    return _angular_average_kernel(grids, Nb, kernel_at_q), Fᴸᴼ
end

"""
Build exactly the enabled microscopic kernel families and apply
[EQ-KSCALE-001](@ref eq-kernel-rescale).  `KernelSet.K` contains normalized
tensors; `KernelSet.qᴷ` contains the exact positive factors restored by every
contraction.

See [Microscopic kernels](@ref theory-kernels) and
[dimensionless kernel scaling](@ref theory-scaling).
"""
function build_kernels(
    p::PhysicalParameters,
    n::NumericalParameters,
    scattering::ScatteringOptions,
    s::ScaleSystem,
    grids::ModelGrids,
    profiles::MaterialProfiles,
    basis::BasisData,
)
    K = Dict{Symbol,Array{ComplexF64,6}}()
    Fᴸᴼ = zeros(ComplexF64, n.N_qz, n.N_b, n.N_b)
    scattering.impurity && (K[:impurity] = impurity_kernel(p, s, grids, profiles, basis))
    scattering.IFR && (K[:IFR] = interface_roughness_kernel(p, s, grids, basis))
    scattering.acoustic && (K[:acoustic] = acoustic_kernel(p, s, grids, basis))
    if scattering.alloy
        p.ΔV_alloy === nothing &&
            throw(ArgumentError("alloy enabled but ΔV_alloy is missing"))
        p.Ω₀ === nothing && throw(ArgumentError("alloy enabled but Ω₀ is missing"))
        K[:alloy] = alloy_kernel(p, s, grids, profiles, basis)
    end
    if scattering.LO
        K[:LO], Fᴸᴼ = lo_phonon_kernel(p, s, grids, basis)
    end
    requested = [name for name in _SCATTERING_MECHANISM_ORDER if haskey(K, name)]
    enabled = Symbol[]
    qᴷ = Dict{Symbol,Float64}()
    normalized = Dict{Symbol,Array{ComplexF64,6}}()
    for name in requested
        q = maximum(abs, K[name])
        isfinite(q) || throw(DomainError(q, "kernel $name is non-finite"))
        if q == 0
            @info "zero microscopic kernel is declared disabled" mechanism=name
            continue
        end
        q > 0 || throw(DomainError(q, "kernel scale must be positive"))
        push!(enabled, name)
        qᴷ[name] = q
        normalized[name] = K[name] ./ q
    end
    return KernelSet(normalized, qᴷ, Fᴸᴼ, enabled)
end

"""
    _static_contraction(K, G, wᵏ)

Literal contraction
`Σ[e,m,a,b]=Σ[m′,c,d] wᵏ[m′]K[m,m′,a,c,d,b]G[e,m′,c,d]`.
No diagonal, local, or momentum-averaged approximation is introduced.
"""
function _static_contraction(
    K::AbstractArray{<:Number,6},
    G::AbstractArray{<:Number,4},
    wᵏ::AbstractVector{<:Real},
    qᴷ::Real = 1.0,
)
    NE, Nk, Nb, Nb2 = size(G)
    Nb == Nb2 || throw(DimensionMismatch("G blocks must be square"))
    size(K) == (Nk, Nk, Nb, Nb, Nb, Nb) ||
        throw(DimensionMismatch("kernel order must be (m,m′,a,c,d,b)"))
    length(wᵏ) == Nk || throw(DimensionMismatch("radial weights differ"))
    Σ = zeros(ComplexF64, NE, Nk, Nb, Nb)
    for e = 1:NE, m = 1:Nk, a = 1:Nb, b = 1:Nb
        value = 0.0 + 0.0im
        for mp = 1:Nk, c = 1:Nb, d = 1:Nb
            value += wᵏ[mp] * K[m, mp, a, c, d, b] * G[e, mp, c, d]
        end
        Σ[e, m, a, b] = qᴷ * value
    end
    return Σ
end

function _lo_contraction(
    K::AbstractArray{<:Number,6},
    green::GreenState,
    problem::NEGFProblem,
    qᴷ::Real,
)
    Gˡ₊ = apply_energy_shift(problem.W₊ᴸᴼ, green.Gˡ)
    Gˡ₋ = apply_energy_shift(problem.W₋ᴸᴼ, green.Gˡ)
    Gᵍ₊ = apply_energy_shift(problem.W₊ᴸᴼ, green.Gᵍ)
    Gᵍ₋ = apply_energy_shift(problem.W₋ᴸᴼ, green.Gᵍ)
    emission_in = _static_contraction(K, Gˡ₊, problem.grids.wᵏ, qᴷ)
    absorption_in = _static_contraction(K, Gˡ₋, problem.grids.wᵏ, qᴷ)
    emission_out = _static_contraction(K, Gᵍ₋, problem.grids.wᵏ, qᴷ)
    absorption_out = _static_contraction(K, Gᵍ₊, problem.grids.wᵏ, qᴷ)
    Nᴸᴼ = _lo_population(problem, green, emission_out, absorption_out)
    Σˡ = (Nᴸᴼ+1) .* emission_in .+ Nᴸᴼ .* absorption_in
    Σᵍ = (Nᴸᴼ+1) .* emission_out .+ Nᴸᴼ .* absorption_out
    return Σˡ, Σᵍ
end

function _symmetric_energy_order(e::Int, NE::Int)
    order = Int[]
    sizehint!(order, NE - 1)
    for δ = 1:max(e-1, NE-e)
        e - δ ≥ 1 && push!(order, e - δ)
        e + δ ≤ NE && push!(order, e + δ)
    end
    return order
end

function _direct_hilbert_ordered(Γ, ε, wᴱ; reverse_order::Bool = false)
    NE, Nk, Nb, _ = size(Γ)
    Λ = zeros(ComplexF64, size(Γ))
    for e = 1:NE
        order = _symmetric_energy_order(e, NE)
        reverse_order && reverse!(order)
        for ep in order
            factor = wᴱ[ep] / (2π * (ε[e] - ε[ep]))
            for m = 1:Nk, a = 1:Nb, b = 1:Nb
                Λ[e, m, a, b] += factor * Γ[ep, m, a, b]
            end
        end
    end
    return Λ
end

"""
    direct_hilbert_transform(Γ, ε, wᴱ; return_roundoff=false)

Direct `O(N_E²)` principal-value quadrature of every matrix element.  The
diagonal energy cell is omitted; the energy grid is never treated as cyclic.
Terms are accumulated in symmetric distance order around the omitted cell.
With `return_roundoff=true`, also return the relative difference obtained by
reversing that fixed order, as required by the Float64 cancellation contract.

Implements [EQ-HILBERT-001](@ref eq-direct-hilbert).
See [Microscopic kernels](@ref theory-kernels) for the causal
principal-value relation.
"""
function direct_hilbert_transform(
    Γ::AbstractArray{<:Number,4},
    ε::AbstractVector{<:Real},
    wᴱ::AbstractVector{<:Real};
    return_roundoff::Bool = false,
)
    NE, Nk, Nb, Nb2 = size(Γ)
    Nb == Nb2 || throw(DimensionMismatch("Γ blocks must be square"))
    length(ε) == NE == length(wᴱ) || throw(DimensionMismatch("energy data differ"))
    Λ = _direct_hilbert_ordered(Γ, ε, wᴱ)
    return_roundoff || return Λ
    Λreverse = _direct_hilbert_ordered(Γ, ε, wᴱ; reverse_order = true)
    r_roundoff = norm(Λ - Λreverse) / (norm(Λ) + 1e-14)
    return Λ, r_roundoff
end

"""
    product_integration_hilbert_matrix(ε)

Analytic product-integration operator for
`PV ∫ Γ(x)/(2π*(E-x)) dx`. Γ is piecewise linear between energy nodes.
Two explicit exterior zero nodes, one adjacent spacing beyond each end,
define a linear edge taper, followed by zero continuation. This boundary
model keeps the transform finite at the first/last sampled nodes; it does
not extrapolate a physical spectral tail. Window convergence must still be
checked away from the taper.

Subtracting Γ(E) analytically cancels both sides of each local singularity.
No diagonal cell is omitted. Works on nonuniform strictly increasing grids.
The returned dense operator can be reused for all columns and SCBA steps.
"""
function product_integration_hilbert_matrix(ε::AbstractVector{<:Real})
    N = length(ε)
    N >= 2 || throw(ArgumentError("at least two energy nodes are required"))
    all(isfinite, ε) && all(>(0), diff(ε)) ||
        throw(ArgumentError("energy nodes must be finite and strictly increasing"))
    x = vcat(ε[1] - (ε[2] - ε[1]), Float64.(ε), ε[end] + (ε[end] - ε[end-1]))
    H = zeros(Float64, N, N)
    for i = 1:N
        E = Float64(ε[i])
        # Integral of the subtracted constant on the entire support.
        H[i, i] += log((E - x[1]) / (x[end] - E))
        for j = 1:(N+1)
            left, right = x[j], x[j+1]
            width = right - left
            # Integral of -(Γ_right-Γ_left)/width, including boundary zeros.
            2 <= j <= N+1 && (H[i, j-1] += 1)
            2 <= j+1 <= N+1 && (H[i, j] -= 1)
            # On an interval touching E the residual numerator is exactly
            # linear in x-E; its logarithmic coefficient vanishes identically.
            (E == left || E == right) && continue
            logarithm = log(abs((E - left) / (E - right)))
            2 <= j <= N+1 && (H[i, j-1] += (right - E) / width * logarithm)
            2 <= j+1 <= N+1 && (H[i, j] += (E - left) / width * logarithm)
            H[i, i] -= logarithm
        end
    end
    return H ./ (2π)
end

"""
    product_integration_hilbert_transform(Γ, ε, wᴱ; return_roundoff=false, operator=nothing)

Apply analytic piecewise-linear PV integration to every tensor column.
`wᴱ` is checked for positive finite values but is not used as an additional
factor: the product-integration matrix already contains the integration
measure. `operator` permits reuse of `product_integration_hilbert_matrix(ε)`.
The optional roundoff diagnostic compares opposite column summation orders;
it measures arithmetic sensitivity, not continuum discretization error.
"""
function product_integration_hilbert_transform(
    Γ::AbstractArray{<:Number,4},
    ε::AbstractVector{<:Real},
    wᴱ::AbstractVector{<:Real};
    return_roundoff::Bool = false,
    operator = nothing,
)
    N, _, Nb, Nb2 = size(Γ)
    Nb == Nb2 || throw(DimensionMismatch("Γ blocks must be square"))
    length(ε) == N == length(wᴱ) || throw(DimensionMismatch("energy data differ"))
    all(isfinite, wᴱ) && all(>(0), wᴱ) || throw(ArgumentError("invalid energy weights"))
    H = operator === nothing ? product_integration_hilbert_matrix(ε) : operator
    size(H) == (N, N) || throw(DimensionMismatch("Hilbert operator differs"))
    columns = reshape(Γ, N, :)
    Λ = reshape(ComplexF64.(H * columns), size(Γ))
    return_roundoff || return Λ
    reverse_result = H[:, end:-1:1] * columns[end:-1:1, :]
    residual = norm(vec(Λ) - vec(reverse_result)) / max(norm(Λ), eps(Float64))
    return Λ, residual
end

"""
    retarded_self_energy(Σˡ, Σᵍ, grids; return_roundoff=false)

Build `Γ=i(Σᵍ-Σˡ)` and `Σᴿ=H[Γ]-iΓ/2` with the direct principal-value
quadrature.  Returns `(Σᴿ,Γ)`. With `return_roundoff=true`, returns
`(Σᴿ,Γ,r_roundoff)`, where the last value compares symmetric and reversed
accumulation orders.

Implements [EQ-SELFENERGY-001](@ref eq-retarded-selfenergy).
See [Microscopic kernels](@ref theory-kernels) and
[the SCBA fixed point](@ref theory-scba).
"""
function retarded_self_energy(
    Σˡ::AbstractArray{<:Number,4},
    Σᵍ::AbstractArray{<:Number,4},
    grids::ModelGrids;
    return_roundoff::Bool = false,
)
    size(Σˡ) == size(Σᵍ) || throw(DimensionMismatch("self-energy components differ"))
    Γ = im .* (Σᵍ .- Σˡ)
    if return_roundoff
        Λ, r_roundoff =
            direct_hilbert_transform(Γ, grids.ε, grids.wᴱ; return_roundoff = true)
        return ComplexF64.(Λ .- 0.5im .* Γ), ComplexF64.(Γ), r_roundoff
    end
    Λ = direct_hilbert_transform(Γ, grids.ε, grids.wᴱ)
    return ComplexF64.(Λ .- 0.5im .* Γ), ComplexF64.(Γ)
end

function _scattering_candidate(
    problem::NEGFProblem,
    green::GreenState;
    return_roundoff::Bool = false,
)
    candidates = Dict{Symbol,SelfEnergyFamily}()
    r_roundoff = 0.0
    for mechanism in problem.kernels.enabled
        if mechanism === :electron_electron
            family, residual = sppa_self_energy(problem, green; return_roundoff = true)
            candidates[mechanism] = family
            r_roundoff = max(r_roundoff, residual)
            continue
        end
        K = problem.kernels.K[mechanism]
        qᴷ = problem.kernels.qᴷ[mechanism]
        if mechanism === :LO
            Σˡ, Σᵍ = _lo_contraction(K, green, problem, qᴷ)
        else
            Σˡ = _static_contraction(K, green.Gˡ, problem.grids.wᵏ, qᴷ)
            Σᵍ = _static_contraction(K, green.Gᵍ, problem.grids.wᵏ, qᴷ)
        end
        if return_roundoff
            Σᴿ, _, rH = retarded_self_energy(Σˡ, Σᵍ, problem.grids; return_roundoff = true)
            r_roundoff = max(r_roundoff, rH)
        else
            Σᴿ, _ = retarded_self_energy(Σˡ, Σᵍ, problem.grids)
        end
        candidates[mechanism] = SelfEnergyFamily(Σᴿ, Σˡ, Σᵍ)
    end
    return return_roundoff ? (candidates, r_roundoff) : candidates
end

function _sum_selfenergies(
    families::AbstractDict{Symbol,SelfEnergyFamily},
    shape::NTuple{4,Int},
)
    total = SelfEnergyFamily(
        zeros(ComplexF64, shape),
        zeros(ComplexF64, shape),
        zeros(ComplexF64, shape),
    )
    for family in values(families)
        total.Σᴿ .+= family.Σᴿ
        total.Σˡ .+= family.Σˡ
        total.Σᵍ .+= family.Σᵍ
    end
    return total
end

function _mix(old::SelfEnergyFamily, candidate::SelfEnergyFamily, α::Real)
    0 < α ≤ 1 || throw(ArgumentError("mixing factor must lie in (0,1]"))
    return SelfEnergyFamily(
        (1 - α) .* old.Σᴿ .+ α .* candidate.Σᴿ,
        (1 - α) .* old.Σˡ .+ α .* candidate.Σˡ,
        (1 - α) .* old.Σᵍ .+ α .* candidate.Σᵍ,
    )
end

function _mix_mechanisms(
    old::Dict{Symbol,SelfEnergyFamily},
    candidate::Dict{Symbol,SelfEnergyFamily},
    α::Real,
)
    return Dict(name => _mix(old[name], candidate[name], α) for name in keys(candidate))
end
