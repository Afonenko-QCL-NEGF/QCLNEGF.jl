const COLLISION_CHANNELS = (:LO,:acoustic,:impurity,:IFR,:alloy)

struct LocalCollisionKernels{LO,IMP}
    coordinate_frame::Matrix{ComplexF64}
    weights_k::Vector{Float64}
    acoustic::Vector{Float64}
    alloy::Vector{Float64}
    interface_probes::Matrix{ComplexF64}
    interface_pair::Matrix{Float64}
    lo_covariance::LO
    impurity_covariance::IMP
    lo_plus::SparseMatrixCSC{Float64,Int}
    lo_minus::SparseMatrixCSC{Float64,Int}
    bose::Float64
end

"""Build local kernels from device-supplied spatial vertices and energy translations.

Covariances use (z,z′,incoming k,outgoing k), before multiplication by wk.
Interface rows are interpolation probes in the physical coordinate period.
"""
function LocalCollisionKernels(q, wk, acoustic, alloy, interface_sites, interface_pair,
        lo_covariance, impurity_covariance, lo_plus, lo_minus, bose)
    nz,r = size(q)
    nk = length(wk)
    length(acoustic) == length(alloy) == nz && size(interface_sites,2)==nz ||
        throw(DimensionMismatch("spatial strengths and probes require the physical coordinate grid"))
    size(interface_pair) == (nk,nk) || throw(DimensionMismatch("interface covariance requires both radial axes"))
    size(lo_covariance) == size(impurity_covariance) == (nz,nz,nk,nk) ||
        throw(DimensionMismatch("spatial covariances require two site and two radial axes"))
    size(lo_plus) == size(lo_minus) && size(lo_plus,1)==size(lo_plus,2) ||
        throw(DimensionMismatch("LO transfer pair requires a square common energy measure"))
    all(>(0),wk) && bose>=0 || throw(ArgumentError("negative physical weights"))
    rank_error = norm(q' * q - I)/sqrt(r)
    rank_error < 1e-10 || throw(ArgumentError("local coordinate frame must be orthonormal"))
    return LocalCollisionKernels(Matrix{ComplexF64}(q), collect(Float64,wk),
        collect(Float64,acoustic),collect(Float64,alloy),Matrix{ComplexF64}(interface_sites*q),
        Matrix{Float64}(interface_pair),lo_covariance,impurity_covariance,
        sparse(Float64.(lo_plus)),sparse(Float64.(lo_minus)),Float64(bose))
end

function contract_spatial!(dest, fields, cov, weights, tile_index, k_out)
    fill!(dest,0)
    @views for kin in eachindex(weights)
        phase = fields[:,:,kin,tile_index]
        kernel = cov[:,:,kin,k_out]
        w = weights[kin]
        @inbounds @simd for i in eachindex(dest)
            dest[i] += w*kernel[i]*phase[i]
        end
    end
    return dest
end

"""Recompute all five local scattering mechanisms in bounded energy tiles.

Destination ordering is rank²×LG×E×k×channel.
The old normalized G and emitted Sigma are packed in the same physical measure.
"""
function scatter_packed!(dest, green, ker::LocalCollisionKernels; tile=8, observe=_->nothing,
        lo_occupation_transfer=nothing)
    q = ker.coordinate_frame
    nz,r = size(q)
    r2,parts,ne,nk = size(green)
    r2 == r^2 && parts==2 && nk==length(ker.weights_k) && size(ker.lo_plus,1)==ne ||
        throw(DimensionMismatch("packed occupation shape differs from geometry or measures"))
    size(dest) == (r2,2,ne,nk,length(COLLISION_CHANNELS)) || throw(DimensionMismatch("wrong output array"))
    tile > 0 || throw(ArgumentError("nonpositive energy tile"))
    tile = min(tile,ne)
    if !isnothing(lo_occupation_transfer)
        length(lo_occupation_transfer)==2 && all(T -> size(T)==(ne,ne) && all(isfinite,T) && all(>=(0),T),
            lo_occupation_transfer) || throw(ArgumentError("custom LO occupation matrices must be finite and nonnegative"))
    end
    site = ker.interface_probes
    faces = size(site,1)
    zbuf = zeros(ComplexF64,nz,nz)
    reduced = zeros(ComplexF64,r,r)
    qtemp = zeros(ComplexF64,nz,r)
    site_temp = zeros(ComplexF64,faces,r)
    elastic = Array{ComplexF64}(undef,nz,nz,nk,tile)
    lo = similar(elastic)
    occupation = zeros(ComplexF64,r,r)
    z_charge = zeros(Float64,nz,tile)
    face_charge = zeros(Float64,faces,nk,tile)
    charge = zeros(Float64,faces)
    lplus,lminus,nbose = ker.lo_plus,ker.lo_minus,ker.bose
    @views for c in 1:2
        incoming,reverse_source = c==1 ? (lplus,lminus) : (lminus,lplus)
        shifted = Array{Float64}(undef,r2,ne,nk)
        for k in 1:nk
            if isnothing(lo_occupation_transfer)
                shifted[:,:,k] .= (1+nbose).*(green[:,c,:,k]*transpose(incoming))
                if nbose>0
                    shifted[:,:,k] .+= nbose.*(green[:,c,:,k]*transpose(reverse_source))
                end
            else
                shifted[:,:,k] .= green[:,c,:,k]*transpose(lo_occupation_transfer[c])
            end
        end
        for start in 1:tile:ne
            span = start:min(start+tile-1,ne)
            fill!(z_charge,0)
            for (i,e) in enumerate(span), kin in 1:nk
                unpack_hermitian!(occupation,green[:,c,e,kin])
                mul!(qtemp,q,occupation)
                mul!(elastic[:,:,kin,i],qtemp,q')
                z_charge[:,i] .+= ker.weights_k[kin].*real.(diag(elastic[:,:,kin,i]))
                mul!(site_temp,site,occupation)
                @inbounds for f in 1:faces
                    face_charge[f,kin,i] = real(dot(site[f,:],site_temp[f,:]))
                end
                unpack_hermitian!(occupation,shifted[:,e,kin])
                mul!(qtemp,q,occupation)
                mul!(lo[:,:,kin,i],qtemp,q')
            end
            for (i,e) in enumerate(span)
                for (channel,strength) in ((2,ker.acoustic),(5,ker.alloy))
                    qtemp .= (strength.*z_charge[:,i]).*q
                    mul!(reduced,q',qtemp)
                    for kout in 1:nk
                        pack_hermitian!(dest[:,c,e,kout,channel],reduced)
                    end
                end
                for kout in 1:nk
                    fill!(charge,0)
                    for kin in 1:nk
                        charge .+= ker.weights_k[kin]*ker.interface_pair[kin,kout].*face_charge[:,kin,i]
                    end
                    site_temp .= charge.*site
                    mul!(reduced,site',site_temp)
                    pack_hermitian!(dest[:,c,e,kout,4],reduced)
                    for (channel,field,cov) in ((1,lo,ker.lo_covariance),(3,elastic,ker.impurity_covariance))
                        contract_spatial!(zbuf,field,cov,ker.weights_k,i,kout)
                        mul!(qtemp,zbuf,q)
                        mul!(reduced,q',qtemp)
                        pack_hermitian!(dest[:,c,e,kout,channel],reduced)
                    end
                end
            end
            observe((occupation=c,last_energy=last(span)))
        end
    end
    return dest
end
