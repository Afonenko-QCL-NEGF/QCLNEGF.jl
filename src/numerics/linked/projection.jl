function unpack_ring(packed,rank)
    ne=size(packed,3)
    l,g=zeros(ComplexF64,ne,1,rank,rank),zeros(ComplexF64,ne,1,rank,rank)
    tmp=zeros(ComplexF64,rank,rank)
    @views for e in 1:ne
        unpack_hermitian!(tmp,packed[:,1,e])
        l[e,1,:,:] .= im.*tmp
        unpack_hermitian!(tmp,packed[:,2,e])
        g[e,1,:,:] .= -im.*tmp
    end
    return SelfEnergyFamily(zeros(ComplexF64,size(l)),l,g)
end

function pack_projected!(out, retarded, fields, k)
    scaled_psd=Inf
    rank,_,_,ne=size(fields)
    @views for e in 1:ne
        if !isnothing(retarded)
            retarded[:,:,e,k] .= fields[:,:,1,e]
        end
        for (component, sign) in ((1,-im),(2,im))
            kinetic=sign.*fields[:,:,component+1,e]
            occup=Hermitian((kinetic+kinetic')/2)
            scaled_psd=min(scaled_psd,eigmin(occup)/max(norm(occup),1e-100))
            pack_hermitian!(out[:,component,e,k],occup)
        end
    end
    return scaled_psd
end

function projected_bond_profile(p,result,k,weight)
    q=p.v.coordinate_frame
    moments=result.moments
    measure=p.we.*weight
    j0=p.current_scale
    cur=zeros(p.nz+1,2)
    for component in 1:2
        integrated=zeros(ComplexF64,p.r,p.r)
        @views for e in 1:p.ne
            integrated .+= measure[e].*moments[:,:,component+1,e]
        end
        coordinate=q*integrated*q'
        for i in 1:p.nz-1
            cur[i+1,component]=2j0*p.hops[i+1,k]*real(coordinate[i+1,i])
        end
    end
    left_raw,right_raw=-sum(result.seams[1,2,:]),sum(result.seams[1,1,:])
    gr_left,gr_right=sum(result.seams[3,2,:]),-sum(result.seams[3,1,:])
    cur[1,:].=j0.*[left_raw,gr_left]
    cur[end,:].=j0.*[right_raw,gr_right]
    return cur
end

"""
Drop-in raw occupation stage for a packed linked local source. The
normalization caller can consume the returned occupation, number, holes
and curr just as in the nodal driver.
"""
function linked_projected_stage!(occupation,p,total,anchor;retarded=nothing,
        order=12,subdivisions=1,grading=true,workers=Threads.nthreads(:default),
        selected_rows=false,projector=project_linked,observe=_->nothing)
    size(occupation)==size(total)==(p.r^2,2,p.ne,p.nk) ||
        throw(DimensionMismatch("source and projected packed occupations differ"))
    size(anchor)==(p.r,p.r,p.nk) || throw(DimensionMismatch("causal reference dimensions differ"))
    ne,nk,r=p.ne,p.nk,p.r
    width=diff(p.en)
    hats=[width[1]; width[1:end-1]+width[2:end]; width[end]]./2
    isapprox(hats,p.we;rtol=5e-11,atol=2e-15) ||
        throw(ArgumentError("linked phase projection uses physical energy hat weights"))
    numbers,holes,spectral_defects,psd_checks=ntuple(_->zeros(nk),4)
    pole_counts,sampling_counts=zeros(Int,nk),zeros(Int,nk)
    currents=zeros(p.nz+1,2,nk)
    seam_lesser=fill(NaN,2,nk)
    timings=zeros(nk)
    workers > 0 || throw(ArgumentError("worker count must be positive"))
    workers=min(workers,Threads.nthreads(:default),nk)
    N._production_worker_ranges!(nk,workers) do _, first_k,last_k
        for k in first_k:last_k
            started=time()
            static=anchor[:,:,k]
            layout=CausalGeometry(p.en,tapered_support(p.en),static,
                p.geo.voltage,maximum(abs,p.geo.offsets))
            source=CausalRing(layout,unpack_ring(view(total,:,:,:,k),r))
            sampler=selected_rows ? observed_sampler(p.geo) : sample_operator!
            radial_weight=p.wk[k]*p.spin_degeneracy/(2π)
            projected=projector(source,p.geo,p.h[:,:,k],p.hops[:,k],radial_weight;
                subdivisions,order,grading,sampler)
            result=projected.result
            psd_checks[k]=pack_projected!(occupation,retarded,result.moments,k)
            numbers[k]=sum(result.occupied)
            holes[k]=sum(result.holes)
            spectral_defects[k]=result.defect
            currents[:,:,k].=projected_bond_profile(p,result,k,radial_weight)
            if hasproperty(result,:seam_lesser)
                seam_lesser[:,k].=result.seam_lesser
            end
            pole_counts[k]=length(projected.poles.poles)
            sampling_counts[k]=length(projected.rule.points)
            timings[k]=time()-started
            observe((;k,charge= numbers[k],spectral_defect=spectral_defects[k],
                poles=pole_counts[k],samples=sampling_counts[k],seconds=timings[k]))
        end
    end
    minimum(psd_checks)>=-1e-9 || error("linked quadrature loses packed positivity")
    number=sum(numbers)
    reconstructed=0.0
    @inbounds for k in 1:nk, e in 1:ne
        reconstructed += p.weights[e,k]*sum(view(occupation,1:r,1,e,k))
    end
    isapprox(reconstructed,number;rtol=1e-10,atol=1e-16) ||
        error("packed projection lost charge")
    curr=dropdims(sum(currents;dims=3);dims=3)
    raw=Dict("raw_sheet_cm2"=>number*p.density_scale,
        "raw_holes_cm2"=>sum(holes)*p.density_scale,
        "raw_lambda"=>number/p.target_dimless,
        "raw_outward_left_right_A_cm2"=>[-curr[1,1],curr[end,1]],
        "raw_profile_A_cm2"=>curr[:,1],"raw_greater_profile_A_cm2"=>curr[:,2],
        "per_ring_raw_sheet_cm2"=>numbers.*p.density_scale,
        "spectral_weight_relative_by_ring"=>spectral_defects,"kinetic_PSD_min_relative"=>psd_checks,
        "pole_count_by_ring"=>pole_counts,"sample_count_by_ring"=>sampling_counts,
        "projection_seconds_by_ring"=>timings)
    return (;occupation,retarded,number,holes=sum(holes),curr,raw,
        seam_lesser=all(isfinite,seam_lesser) ? seam_lesser : nothing)
end
