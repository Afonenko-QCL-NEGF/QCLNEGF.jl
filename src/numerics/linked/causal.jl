struct CausalRows
    query::Matrix{Float64}
    kinetic::Matrix{Float64}
    pv::Matrix{Float64}
end

function stable_local_log(edge,z,width)
    edge == 0 && return 0.0
    ratio = width/z
    return abs(ratio) <= 0.5 ? log1p(ratio) : log(abs(edge/z))
end

function nodal_hat_pv(q,a,b,c)
    dl,dr = b-a,c-b
    z = q-b
    iszero(z) && return log(dl/dr)/(2π)
    left,right = q-a,q-c
    return (left/dl*stable_local_log(left,z,dl)+right/dr*stable_local_log(right,z,-dr))/(2π)
end

function CausalRows(energies, shifted)
    nq,nshifts = size(shifted)
    ne = length(energies)
    support=[2first(energies)-energies[2];energies;2last(energies)-energies[end-1]]
    all(diff(support) .> 0) || throw(ArgumentError("non-increasing energies"))
    snapped = Array{Float64}(shifted)
    gap = 64eps(Float64)*max(maximum(abs,energies),1)
    @inbounds for i in eachindex(snapped)
        q = snapped[i]
        high = min(searchsortedfirst(energies,q),ne)
        low = max(high-1,1)
        closest = abs(q-energies[high]) <= abs(q-energies[low]) ? high : low
        if abs(q-energies[closest]) <= gap
            snapped[i] = energies[closest]
        end
    end
    pv = Matrix{Float64}(undef,nq*nshifts,ne)
    @inbounds for e in eachindex(energies), i in eachindex(snapped)
        pv[i,e] = nodal_hat_pv(snapped[i],support[e],support[e+1],support[e+2])
    end
    kinetic = zeros(Float64, nq*nshifts,ne)
    @inbounds for i in eachindex(snapped)
        q = snapped[i]
        if first(support) <= q <= last(support)
            lo = clamp(searchsortedlast(support,q),1,ne+1)
            x = (q-support[lo])/(support[lo+1]-support[lo])
            if lo > 1
                kinetic[i,lo-1] = 1-x
            end
            if lo <= ne
                kinetic[i,lo] = x
            end
        end
    end
    return CausalRows(snapped,kinetic,pv)
end

function sample_packed_ring(plan::CausalRows, occupation, reference)
    r2,parts,ne = size(occupation)
    rank = isqrt(r2)
    rank^2==r2 && parts==2 && size(reference)==(rank,rank) && ne==size(plan.kinetic,2) ||
        throw(DimensionMismatch("packed ring requires rank²×LG×E"))
    return sample_real_ring(plan,transpose(view(occupation,:,1,:)),transpose(view(occupation,:,2,:)),reference)
end

function sample_real_ring(plan,lk,gk,reference)
    rank = size(reference,1)
    l,g = plan.kinetic*lk,plan.kinetic*gk
    dispersive = plan.pv*(lk+gk)
    result = zeros(ComplexF64,rank,rank,3,size(plan.query)...)
    @views for (i,key) in enumerate(CartesianIndices(plan.query))
        e,p = Tuple(key)
        r,d,f = result[:,:,1,e,p],result[:,:,2,e,p],result[:,:,3,e,p]
        unpack_hermitian!(r,dispersive[i,:])
        unpack_hermitian!(d,l[i,:])
        unpack_hermitian!(f,g[i,:])
        d .*= im
        f .*= -im
        r .+= reference .+ (f .- d)./2
    end
    return result
end

function tapered_support(ε)
    length(ε) >= 2 && all(isfinite, ε) && all(>(0), diff(ε)) ||
        throw(ArgumentError("energy support must be finite and strictly increasing"))
    return [2ε[1]-ε[2]; ε; 2ε[end]-ε[end-1]]
end

# PV integral of one triangular hat against 1/(2π(E-x)). The two endpoint
# terms with zero coefficients are assigned their finite limits explicitly.
function hat_pv(E, left, center, right)
    a, b, c = E-left, E-center, E-right
    dl, dr = center-left, right-center
    iszero(b) && return log(dl/dr)/(2π)
    ql, qr = dl/b, -dr/b
    hl = iszero(a) ? 0.0 : a/dl * (abs(ql) <= 0.5 ? log1p(ql) : log(abs(a/b)))
    hr = iszero(c) ? 0.0 : c/dr * (abs(qr) <= 0.5 ? log1p(qr) : log(abs(c/b)))
    return (hl + hr)/(2π)
end

function point_sources(energy, point)
    first(energy) <= point <= last(energy) || return (0, 0, 0.0)
    left = searchsortedlast(energy, point)
    left == lastindex(energy) && return (left, left, 0.0)
    right = left+1
    return (left, right, (point-energy[left])/(energy[right]-energy[left]))
end

@inline function hermitian_width(family, e, m, a, b)
    difference = family.Σᵍ[e, m, a, b] - family.Σˡ[e, m, a, b]
    transpose = family.Σᵍ[e, m, b, a] - family.Σˡ[e, m, b, a]
    return 0.5im*(difference-conj(transpose))
end

"""Geometry for real-axis PV evaluation over the finite chain."""
struct CausalGeometry
    energy::Vector{Float64}
    support::Vector{Float64}
    constant::Matrix{ComplexF64}
    drop::Float64
    periods::Int
end

"""Pack the upper triangular local width of one ring for small real BLAS batches."""
struct CausalRing
    geometry::CausalGeometry
    local_sum::SelfEnergyFamily
    gamma::Matrix{Float64}
    triangle::Vector{Tuple{Int,Int}}
end

function CausalRing(geometry, local_sum)
    ne = length(geometry.energy)
    nb = size(geometry.constant, 1)
    triangle = [(a,b) for b in 1:nb for a in 1:b]
    gamma = Matrix{Float64}(undef, ne, 2length(triangle))
    @inbounds for (i,(a,b)) in enumerate(triangle), e in 1:ne
        gamma[e,2i-1], gamma[e,2i] = reim(hermitian_width(local_sum,e,1,a,b))
    end
    return CausalRing(geometry, local_sum, gamma, triangle)
end

function kinetic_row(g, point)
    energy, support = g.energy, g.support
    first(support) < point < first(energy) &&
        return (0, 1, (point-first(support))/(first(energy)-first(support)))
    last(energy) < point < last(support) &&
        return (length(energy), 0, (point-last(energy))/(last(support)-last(energy)))
    return point_sources(energy, point)
end

@inline function sources_at(field, row, k, a, b)
    left, right, fraction = row
    low = iszero(left) ? 0.0im : field[left,k,a,b]
    iszero(fraction) && return low
    high = iszero(right) ? 0.0im : field[right,k,a,b]
    return (1-fraction)*low + fraction*high
end

function source_triplets(ring, points; block_size=192)
    g = ring.geometry
    ns, nb, ne = length(points), size(g.constant, 1), length(g.energy)
    local_field = N._empty_selfenergy(ns, 1, nb)
    for first_q in 1:block_size:ns
        last_q = min(first_q+block_size-1, ns)
        indices = first_q:last_q
        H = Matrix{Float64}(undef, length(indices), ne)
        @inbounds for e in 1:ne, (row, q) in enumerate(indices)
            H[row,e] = hat_pv(points[q], g.support[e], g.support[e+1], g.support[e+2])
        end
        result = H * ring.gamma
        fill_triplets!(local_field, ring, points, indices, result)
    end
    return local_field
end

function fill_triplets!(local_field, ring, points, indices, pv)
    g, source = ring.geometry, ring.local_sum
    @inbounds for (row,q) in enumerate(indices)
        kinetic = kinetic_row(g, points[q])
        for (i,(a,b)) in enumerate(ring.triangle)
            rab = g.constant[a,b]+complex(pv[row,2i-1],pv[row,2i])
            lab = sources_at(source.Σˡ,kinetic,1,a,b)
            gba = sources_at(source.Σᵍ,kinetic,1,b,a)
            gab = sources_at(source.Σᵍ,kinetic,1,a,b)
            lba = sources_at(source.Σˡ,kinetic,1,b,a)
            anti = (gab-lab-conj(gba-lba))/4
            local_field.Σᴿ[q,1,a,b] = rab+anti
            local_field.Σᴿ[q,1,b,a] = conj(rab)-conj(anti)
            local_field.Σˡ[q,1,a,b] = lab
            local_field.Σˡ[q,1,b,a] = lba
            local_field.Σᵍ[q,1,a,b] = gab
            local_field.Σᵍ[q,1,b,a] = gba
        end
    end
    return local_field
end

function shifted_rows(g, x)
    n = length(x)
    offsets = -g.periods:g.periods
    translated = [point + shift*g.drop for point in x, shift in offsets]
    # Equal translated energies share the same PV evaluation.
    points = unique(vec(translated))
    by_value = Dict(value => i for (i, value) in enumerate(points))
    rows = [[(by_value[translated[e,j]]=>1.0,) for e in 1:n] for j in axes(translated,2)]
    return points, rows
end
