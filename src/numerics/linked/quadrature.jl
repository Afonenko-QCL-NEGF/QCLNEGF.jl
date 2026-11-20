# Positive Gauss–Legendre nodes on the reference interval.
function gauss_rule(order)
    order >= 2 || throw(ArgumentError("use Gauss order at least two"))
    steps = collect(1:order-1)
    gauss = eigen(SymTridiagonal(zeros(order), steps ./ sqrt.(4 .* steps.^2 .- 1)))
    return (;nodes=gauss.values, weights=2 .* gauss.vectors[1, :].^2)
end

function inverse_samples(ring, tail, h, energy)
    points,rows=shifted_rows(ring.geometry,energy)
    source=source_triplets(ring,points)
    modes=size(h,1)
    rank=size(first(tail.bases),1)
    temp=zeros(ComplexF64,rank,modes)
    result=Vector{Matrix{ComplexF64}}(undef,length(energy))
    for i in eachindex(energy)
        D=energy[i].*tail.metric-h
        for (j,basis) in enumerate(tail.bases)
            q=first(rows[j][i]).first
            @views mul!(temp,source.Σᴿ[q,1,:,:],basis)
            mul!(D,adjoint(basis),temp,-1,1)
        end
        result[i]=D
    end
    return result
end

struct Crossing
    eigindex::Int
    lo::Float64
    hi::Float64
    f_lo::Float64
    f_hi::Float64
end

hermitian_part(A) = Hermitian((A+A')/2)

function translated_knots(ring, extra_splits)
    ε=ring.geometry.energy
    lo,hi=extrema(ε)
    translated=sort!([y for j in -ring.geometry.periods:ring.geometry.periods
        for x in ring.geometry.support for y in (x-j*ring.geometry.drop,) if lo<=y<=hi])
    nodes=Float64[]
    threshold=64eps(Float64)*max(abs(lo),abs(hi),1)
    for e in [lo; translated; hi]
        (isempty(nodes) || e-last(nodes)>threshold) && push!(nodes,e)
    end
    last(nodes)<hi && (nodes[end]=hi)
    expanded=vec([[lo+(hi-lo)*q/extra_splits for q in 0:extra_splits-1] for
        (lo,hi) in zip(nodes[1:end-1],nodes[2:end])])
    return [vcat(expanded...); nodes[end]]
end

function crossings(inverse, points)
    eigs=map(A->eigvals(hermitian_part(A)),inverse)
    found=Crossing[]
    for i in 1:length(points)-1, n in eachindex(first(eigs))
        lo,hi=eigs[i][n],eigs[i+1][n]
        if signbit(lo)!=signbit(hi)
            push!(found,Crossing(n,points[i],points[i+1],lo,hi))
        end
    end
    return found
end

"""
Locate narrow-pole candidates by zero crossings in each ordered eigenvalue of
the Hermitian inverse Dyson pencil on a partition of shifted kinetic knots.
The contour and independent subdivision rule test the resulting pole list.
"""
function scan_poles(ring, tail, h; subdivisions=2, find_steps=28)
    points=translated_knots(ring,subdivisions)
    found=crossings(inverse_samples(ring,tail,h,points),points)
    for _ in 1:find_steps
        mid=[(r.lo+r.hi)/2 for r in found]
        matrices=inverse_samples(ring,tail,h,mid)
        found=map(zip(found,mid,matrices)) do (r,e,D)
            f=eigvals(hermitian_part(D))[r.eigindex]
            if signbit(f)==signbit(r.f_lo)
                Crossing(r.eigindex,e,r.hi,f,r.f_hi)
            else
                Crossing(r.eigindex,r.lo,e,r.f_lo,f)
            end
        end
    end
    grid=ring.geometry.energy
    radii=map(found) do r
        e=(r.f_hi*r.lo-r.f_lo*r.hi)/(r.f_hi-r.f_lo)
        isfinite(e) ? e : (r.lo+r.hi)/2
    end
    widths=map(zip(radii,found,inverse_samples(ring,tail,h,radii))) do (e,r,D)
        λ=eigen(hermitian_part(D))
        v=view(λ.vectors,:,r.eigindex)
        derivative=(r.f_hi-r.f_lo)/(r.hi-r.lo)
        loss=abs(imag(dot(v,D*v)))
        width=loss/max(abs(derivative),1e-15)
        (;energy=e,width,loss,slope=derivative,
            bracket_width=r.hi-r.lo,residual=abs(λ.values[r.eigindex]))
    end
    return (;points,poles=sort!(widths;by=q->q.energy))
end

function integrate_panel(lo,hi,pole,rule; grading=true)
    x,w=rule.nodes,rule.weights
    y=(x.+1)./2
    shape=grading ? y.^2 .* (3 .- 2y) : y
    derivative=grading ? 6y.*(1 .- y) : ones(length(y))
    if isnothing(pole)
        e=lo .+ (hi-lo).*shape
        weight=(hi-lo).*w./2 .* derivative
    else
        center=pole.energy
        width=max(pole.width,8eps(Float64)*max(abs(center),1))
        α,β=asinh((lo-center)/width),asinh((hi-center)/width)
        p=α .+ (β-α).*shape
        e=center .+ width.*sinh.(p)
        weight=width*(β-α).*cosh.(p).*w./2 .* derivative
    end
    # The transformed rule's mass on each disjoint panel is exact.
    weight .*= (hi-lo)/sum(weight)
    return e,weight
end

function resonance_rule(ring, scan, order; grading=true)
    ε=ring.geometry.energy
    centers=[p.energy for p in scan.poles]
    bisectors=[(lo+hi)/2 for (lo,hi) in zip(centers[1:end-1],centers[2:end])]
    refined=sort!(unique!([scan.points; centers; bisectors]))
    points,weights=Float64[],Float64[]
    cell=Int[]
    rule=gauss_rule(order)
    for (lo,hi) in zip(refined[1:end-1],refined[2:end])
        mid=(lo+hi)/2
        nearby=filter(p->lo<=p.energy<=hi,scan.poles)
        # A narrow pole strictly inside a shifted-knot interval is itself a
        # panel boundary. Apply the sinh map on each neighboring half.
        pole=isempty(nearby) ? nothing : nearby[argmin([abs(p.energy-mid) for p in nearby])]
        e,w=integrate_panel(lo,hi,pole,rule;grading)
        append!(points,e)
        append!(weights,w)
        append!(cell,fill(clamp(searchsortedlast(ε,mid),1,length(ε)-1),length(e)))
    end
    abs(sum(weights)-(last(ε)-first(ε)))<=32eps(Float64)*(last(ε)-first(ε)) ||
        error("positive quadrature changed energy measure")
    minimum(weights)>0 || error("adaptive quadrature produced a nonpositive panel weight")
    return (;points,weights,cell)
end

function integrate_linked(ring,tail,h,hops,rule,wk; chunksize=256, project=false,
        sampler=sample_operator!)
    ε=ring.geometry.energy
    nc=length(ε)-1
    rank=size(first(tail.bases),1)
    ws=TailWorkspace(tail)
    spectrum=zeros(ComplexF64,rank,rank)
    occupied,holes=zeros(nc),zeros(nc)
    seams=zeros(3,2,nc)
    seam_lesser=zeros(2)
    moments=project ? zeros(ComplexF64,rank,rank,3,nc+1) : nothing
    widths=diff(ε)
    energy_hats=[widths[1]; widths[1:end-1].+widths[2:end]; widths[end]]./2
    for first_i in 1:chunksize:length(rule.points)
        indices=first_i:min(first_i+chunksize-1,length(rule.points))
        x=rule.points[indices]
        p,rows=shifted_rows(ring.geometry,x)
        source=source_triplets(ring,p)
        results=sampler(ws,tail,x,h,hops,source,rows)
        for (q,r) in enumerate(indices)
            w,c=rule.weights[r],rule.cell[r]
            gr=view(results.gr,q,:,:)
            spectrum .+= (im*w/(2π)).*(gr-gr')
            occupied[c]+=w*wk*imag(tr(view(results.gl,q,:,:)))
            holes[c]-=w*wk*imag(tr(view(results.gg,q,:,:)))
            @views seams[:,:,c] .+= w*wk.*results.seams[:,:,q]
            @views seam_lesser .+= w.*results.seam_lesser[:,q]
            if project
                t=(x[q]-ε[c])/widths[c]
                for (tag,field) in enumerate((results.gr,results.gl,results.gg))
                    @views begin
                        moments[:,:,tag,c].+=(w*(1-t)/energy_hats[c]).*field[q,:,:]
                        moments[:,:,tag,c+1].+=(w*t/energy_hats[c+1]).*field[q,:,:]
                    end
                end
            end
        end
    end
    metric=spectral_completeness(tail)
    if project
        restored=sum(energy_hats[e]*wk*imag(tr(view(moments,:,:,2,e))) for e in eachindex(ε))
        isapprox(restored,sum(occupied);rtol=1e-11,atol=1e-16) ||
            error("linked nodal mass projection changed the resolved sheet charge")
    end
    return (;spectrum,occupied,holes,seams,seam_lesser,moments,defect=norm(spectrum-metric)/norm(metric))
end

"""Resolve all three raw linked Green tensors and project them with positive affine hats."""
function project_linked(ring,tail,h,hops,wk;subdivisions=1,order=12,grading=true,
        sampler=sample_operator!)
    poles=scan_poles(ring,tail,h;subdivisions)
    rule=resonance_rule(ring,poles,order;grading)
    result=integrate_linked(ring,tail,h,hops,rule,wk;sampler,project=true)
    return (;result,poles,rule)
end

function economical_rule(ring, scan, orders, cutoffs; grading=true)
    energy = ring.geometry.energy
    poles = scan.poles
    centers = [pole.energy for pole in poles]
    midpoints = [(l+r)/2 for (l,r) in zip(centers[1:end-1],centers[2:end])]
    cut = sort!(unique!([scan.points; centers; midpoints]))
    rules = Dict(order=>gauss_rule(order) for order in orders)
    points,weights,cell = Float64[],Float64[],Int[]
    counts = Dict(order=>0 for order in keys(rules))
    for (lo,hi) in zip(cut[1:end-1],cut[2:end])
        width,mid = hi-lo,(lo+hi)/2
        near = filter(pole->lo<=pole.energy<=hi,poles)
        pole = isempty(near) ? nothing : near[argmin([abs(pole.energy-mid) for pole in near])]
        distance = minimum(max(lo-e,e-hi,0) for e in centers;init=Inf)
        order = distance>=cutoffs[1]*width ? orders[1] : (distance>=cutoffs[2]*width ? orders[2] : orders[3])
        e,w = integrate_panel(lo,hi,pole,rules[order];grading)
        append!(points,e)
        append!(weights,w)
        append!(cell,fill(clamp(searchsortedlast(energy,mid),1,length(energy)-1),order))
        counts[order] += 1
    end
    isapprox(sum(weights),last(energy)-first(energy);rtol=32eps(Float64),atol=0) || error("measure mismatch")
    minimum(weights)>0 || error("nonpositive weight")
    return (;points,weights,cell,counts)
end
"""Project with the validated distance-tiered GL9/13/18 rule on two knot subdivisions."""
function project_linked_economical(ring,tail,h,hops,wk;subdivisions=2,order=18,grading=true,
        sampler=sample_operator!)
    order == 18 || throw(ArgumentError("distance-tiered linked projection uses order 18 near poles"))
    poles=scan_poles(ring,tail,h;subdivisions)
    rule=economical_rule(ring,poles,(9,13,18),(5,2);grading)
    result=integrate_linked(ring,tail,h,hops,rule,wk;sampler,project=true)
    return (;result,poles,rule)
end
