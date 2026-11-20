"""Static shape and coordinate observation for one finite translated mode chain."""
struct TailGeometry
    modes::Matrix{ComplexF64}
    bases::Vector{Matrix{ComplexF64}}
    metric::Matrix{ComplexF64}
    center_metric::Matrix{ComplexF64}
    monitor::UnitRange{Int}
    center_modes::UnitRange{Int}
    block_ranges::Vector{UnitRange{Int}}
    offsets::Vector{Int}
    voltage::Float64
    reach::Int
end

function TailGeometry(modes, q, block_width, offsets, voltage, reach)
    n_coord, rank = size(q)
    periods = length(offsets)
    isodd(periods) && periods >= 3 || throw(ArgumentError("expected an odd open chain with two wings"))
    block_width > 0 && rank >= block_width && reach >= 1 || throw(ArgumentError("invalid linked dimensions"))
    offsets == collect(-div(periods, 2):div(periods, 2)) || throw(ArgumentError("periods must be centered unit offsets"))
    size(modes) == (periods*n_coord, periods*block_width) || throw(DimensionMismatch("invalid mode frame"))
    u = Matrix{ComplexF64}(modes)
    coordinate_rows = [(j-1)*n_coord+1:j*n_coord for j in 1:periods]
    mode_columns = [(j-1)*block_width+1:j*block_width for j in 1:periods]
    center = div(periods, 2)+1
    center_z = coordinate_rows[center]
    local_phys = u[center_z,:]
    metric = u' * u
    isposdef(Hermitian(metric)) || throw(ArgumentError("linearly dependent modes"))
    return TailGeometry(u, [q' * u[indices,:] for indices in coordinate_rows], metric,
        local_phys' * local_phys, first(center_z)-1:last(center_z)+1, mode_columns[center],
        mode_columns, collect(Int, offsets), Float64(voltage), Int(reach))
end

function centered_basis(geo::TailGeometry)
    return geo.bases[div(length(geo.bases),2)+1]
end

function spectral_completeness(geo::TailGeometry)
    b = centered_basis(geo)
    return b * (geo.metric \ b')
end

shifted_queries(geo::TailGeometry, en) = [x + p*geo.voltage for x in en, p in geo.offsets]

"""One reusable worker. R/L/G are the last axes of the matrix fields."""
struct TailWorkspace
    mode_source::Array{ComplexF64,3}
    d::Matrix{ComplexF64}
    factor::Matrix{ComplexF64}
    mode_green::Array{ComplexF64,3}
    local_green::Array{ComplexF64,3}
    vertex_scratch::Matrix{ComplexF64}
    kinetic_scratch::Matrix{ComplexF64}
    observer_scratch::Matrix{ComplexF64}
    monitor_scratch::Matrix{ComplexF64}
    bond_profile::Matrix{Float64}
    seam_lesser::Vector{Float64}
    triplets::Array{ComplexF64,4}
end

function TailWorkspace(geo::TailGeometry)
    n_mode = size(geo.modes,2)
    rank = size(first(geo.bases),1)
    nm = length(geo.monitor)
    return TailWorkspace(zeros(ComplexF64,n_mode,n_mode,3), zeros(ComplexF64,n_mode,n_mode),
        zeros(ComplexF64,n_mode,n_mode), zeros(ComplexF64,n_mode,n_mode,3),
        zeros(ComplexF64,rank,rank,3), zeros(ComplexF64,rank,n_mode),
        zeros(ComplexF64,n_mode,n_mode), zeros(ComplexF64,rank,n_mode),
        zeros(ComplexF64,nm,n_mode), zeros(Float64,nm-1,2), zeros(2),
        zeros(ComplexF64,rank,rank,3,length(geo.bases)))
end

function assemble_sources!(ws::TailWorkspace, geo::TailGeometry, local_triplets)
    size(local_triplets) == size(ws.triplets) || throw(DimensionMismatch("expected rank×rank×3×periods"))
    fill!(ws.mode_source,0)
    @views for c in axes(local_triplets,3), j in eachindex(geo.bases)
        b = geo.bases[j]
        mul!(ws.vertex_scratch, local_triplets[:,:,c,j], b)
        mul!(ws.mode_source[:,:,c], b', ws.vertex_scratch, 1, 1)
    end
    return ws.mode_source
end

function invert_chain!(out, work, original, geo::TailGeometry)
    copyto!(work,original)
    fill!(out,0)
    @inbounds for i in axes(out,1)
        out[i,i] = 1
    end
    blocks = geo.block_ranges
    @views for (i,g) in enumerate(blocks)
        pivot = lu(work[g,g]; check=false)
        if !issuccess(pivot)
            # A singular block pivot can occur while the whole finite D is invertible.
            # In that case use partial pivoting on the full operator.
            copyto!(out, inv(original))
            return out
        end
        ldiv!(pivot,out[g,:])
        following = i+1:min(length(blocks),i+geo.reach)
        for j in following
            ldiv!(pivot, work[g,blocks[j]])
        end
        for j in following
            gj = blocks[j]
            mul!(out[gj,:], work[gj,g], out[g,:], -1, 1)
            for k in following
                gk = blocks[k]
                mul!(work[gj,gk], work[gj,g], work[g,gk], -1, 1)
            end
        end
    end
    @views for i in reverse(eachindex(blocks)), j in i+1:min(length(blocks),i+geo.reach)
        mul!(out[blocks[i],:],work[blocks[i],blocks[j]],out[blocks[j],:],-1,1)
    end
    return out
end

function observe_profile!(ws::TailWorkspace, geo::TailGeometry, hopping)
    length(hopping) == size(ws.bond_profile,1) || throw(DimensionMismatch("wrong physical faces"))
    @views coordinates = geo.modes[geo.monitor,:]
    @views for c in 1:2
        mul!(ws.monitor_scratch, coordinates, ws.mode_green[:,:,c+1])
        @inbounds for z in eachindex(hopping)
            ws.bond_profile[z,c] = 2*hopping[z]*real(dot(coordinates[z,:],ws.monitor_scratch[z+1,:]))
        end
        if c == 1
            ws.seam_lesser[1] = real(dot(coordinates[1,:], ws.monitor_scratch[2,:]))
            ws.seam_lesser[2] = real(dot(coordinates[end-1,:], ws.monitor_scratch[end,:]))
        end
    end
    return ws.bond_profile
end

"""At one physical energy, assemble all shifted sources and solve the linked chain."""
function sample_node!(ws::TailWorkspace, geo::TailGeometry, energy, h, hopping, local_triplets)
    assemble_sources!(ws,geo,local_triplets)
    @views begin
        ws.d .= energy .* geo.metric .- h .- ws.mode_source[:,:,1]
        r = invert_chain!(ws.mode_green[:,:,1], ws.factor, ws.d, geo)
        for c in 2:3
            mul!(ws.kinetic_scratch,r,ws.mode_source[:,:,c])
            mul!(ws.mode_green[:,:,c],ws.kinetic_scratch,r')
        end
        b = centered_basis(geo)
        for c in 1:3
            mul!(ws.observer_scratch,b,ws.mode_green[:,:,c])
            mul!(ws.local_green[:,:,c],ws.observer_scratch,b')
        end
    end
    observe_profile!(ws,geo,hopping)
    return ws
end

"""Left and right boundary corrections with explicit mixed lesser/greater sources."""
function boundary_parts(ws::TailWorkspace, geo::TailGeometry)
    c = geo.center_modes
    l, r = 1:first(c)-1, last(c)+1:size(ws.d,1)
    outside = [l;r]
    d, source = ws.d, ws.mode_source
    t = -transpose(transpose(d[outside,outside]) \ transpose(d[c,outside]))
    lt, rt = t[:,1:length(l)], t[:,length(l)+1:end]
    nb = length(c)
    parts = zeros(ComplexF64,nb,nb,3,3)
    for (side,(indices,trans)) in enumerate(((l,lt),(r,rt)))
        parts[:,:,1,side] .= -trans * d[indices,c]
        for p in 2:3
            parts[:,:,p,side] .= trans*source[indices,indices,p]*trans' +
                trans*source[indices,c,p] + source[c,indices,p]*trans'
        end
    end
    for p in 2:3
        parts[:,:,p,3] .= lt*source[l,r,p]*rt' + rt*source[r,l,p]*lt'
    end
    return parts
end

"""Julia uniform-quadrature callback over exact shifted source row plans.

Returned seam columns are the raw value, derivative in occupied fraction a,
and derivative in added empty fraction c, with array side order +z, −z.
The separate lesser coherences carry left and right cuts, both along +z,
before hopping factors and quadrature measures.
"""
function sample_operator!(ws::TailWorkspace, geo::TailGeometry, energy, h, hopping, source, rows)
    n = length(energy)
    rank = size(ws.local_green,1)
    fields = (source.Σᴿ,source.Σˡ,source.Σᵍ)
    result = ntuple(_ -> Array{ComplexF64}(undef,n,rank,rank),3)
    seams = zeros(Float64,3,2,n)
    seam_lesser = zeros(2,n)
    length(rows) == length(geo.offsets) || throw(DimensionMismatch("shifted periods do not match geometry"))
    @views for e in eachindex(energy)
        fill!(ws.triplets,0)
        for (shift,plan) in enumerate(rows), (row,weight) in plan[e], comp in 1:3
            ws.triplets[:,:,comp,shift] .+= weight .* fields[comp][row,1,:,:]
        end
        sample_node!(ws,geo,energy[e],h,hopping,ws.triplets)
        for comp in 1:3
            result[comp][e,:,:] .= ws.local_green[:,:,comp]
        end
        raw_p, raw_m = ws.bond_profile[end,1], -ws.bond_profile[1,1]
        seams[:,1,e] .= (raw_p,raw_p,-ws.bond_profile[end,2])
        seams[:,2,e] .= (raw_m,raw_m,ws.bond_profile[1,2])
        seam_lesser[:,e] .= ws.seam_lesser
    end
    return (; gr=result[1],gl=result[2],gg=result[3],seams,seam_lesser)
end

"""Exactly n² real scalars encode a Hermitian matrix; all diagonals come first."""
function pack_hermitian!(dest, mat)
    n = size(mat,1)
    size(mat,2) == n && length(dest) == n*n || throw(DimensionMismatch("Hermitian pack shape"))
    @inbounds for i in 1:n
        dest[i] = real(mat[i,i])
    end
    p = n+1
    @inbounds for j in 2:n, i in 1:j-1
        dest[p],dest[p+1] = reim(mat[i,j])
        p += 2
    end
    return dest
end

function unpack_hermitian!(mat, src)
    n = size(mat,1)
    size(mat,2) == n && length(src) == n*n || throw(DimensionMismatch("Hermitian unpack shape"))
    @inbounds for i in 1:n
        mat[i,i] = src[i]
    end
    p = n+1
    @inbounds for j in 2:n, i in 1:j-1
        z = complex(src[p],src[p+1])
        mat[i,j],mat[j,i] = z,conj(z)
        p += 2
    end
    return mat
end
