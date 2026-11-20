"""
    LinkedProblem(energy, energy_weights, radial_weights, modes, block_width,
                  offsets, voltage, reach, hamiltonian, hopping, anchor, kernels;
                  spin_degeneracy, current_scale, density_scale, target_dimless)

Physical inputs for a finite Q-linked chain in the dimensionless units of the
packed solver. The local coordinate frame and translated LO weights live in
`kernels`. Hamiltonians are mode × mode × k; the hopping array has one entry
per physical cut (including both seams) and radial ring. The causal anchor is
rank × rank × k × channel in `COLLISION_CHANNELS` order.
"""
struct LinkedProblem{K<:LocalCollisionKernels}
    en::Vector{Float64}
    we::Vector{Float64}
    wk::Vector{Float64}
    weights::Matrix{Float64}
    geo::TailGeometry
    h::Array{ComplexF64,3}
    hops::Matrix{Float64}
    const_by_channel::Array{ComplexF64,4}
    v::K
    ne::Int
    nk::Int
    r::Int
    nb::Int
    nz::Int
    spin_degeneracy::Float64
    current_scale::Float64
    density_scale::Float64
    target_dimless::Float64
    shift_error::Float64
end

function LinkedProblem(energy, energy_weights, radial_weights, modes, block_width,
        offsets, voltage, reach, hamiltonian, hopping, anchor, kernels::LocalCollisionKernels;
        spin_degeneracy, current_scale, density_scale, target_dimless)
    tapered_support(energy)
    ne,nk = length(energy),length(radial_weights)
    length(energy_weights)==ne && !isempty(radial_weights) &&
        all(w->isfinite(w) && w>0,energy_weights) && all(w->isfinite(w) && w>0,radial_weights) ||
        throw(ArgumentError("physical weights must be positive and match the axes"))
    radial_weights == kernels.weights_k ||
        throw(ArgumentError("collision and transport radial weights differ"))
    all(s->isfinite(s) && s>0,(spin_degeneracy,current_scale,density_scale,target_dimless)) ||
        throw(ArgumentError("invalid physical scales"))
    isfinite(voltage) || throw(ArgumentError("nonfinite period voltage"))
    q = kernels.coordinate_frame
    nz,r = size(q)
    m = size(modes,2)
    size(hamiltonian)==(m,m,nk) && size(hopping)==(nz+1,nk) &&
        size(anchor)==(r,r,nk,length(COLLISION_CHANNELS)) && size(kernels.lo_plus)==(ne,ne) ||
        throw(DimensionMismatch("linked inputs differ in geometry, channels, or measure"))
    all(isfinite,hopping) || throw(ArgumentError("nonfinite physical hopping"))
    all(isfinite,anchor) || throw(ArgumentError("nonfinite causal anchor"))
    geo = TailGeometry(modes,q,block_width,offsets,voltage,reach)
    @views for k in 1:nk
        block = hamiltonian[:,:,k]
        norm(block-block') <= 1e-12*max(norm(block),1e-300) ||
            throw(ArgumentError("non-Hermitian mode Hamiltonian at ring $k"))
    end
    @views for i in eachindex(geo.block_ranges), j in eachindex(geo.block_ranges)
        abs(i-j)<=geo.reach && continue
        rows,cols=geo.block_ranges[i],geo.block_ranges[j]
        norm(geo.metric[rows,cols])+norm(hamiltonian[rows,cols,:]) <= 1e-12 ||
            throw(ArgumentError("mode geometry exceeds declared half-bandwidth at block ($i,$j)"))
        any(b->norm(b[:,rows])*norm(b[:,cols])>1e-12,geo.bases) &&
            throw(ArgumentError("local sources exceed declared half-bandwidth at block ($i,$j)"))
    end
    plus,minus = kernels.lo_plus,kernels.lo_minus
    forward = energy_weights.*plus
    shift_error = norm(forward-(energy_weights.*minus)')/max(norm(forward),floatmin(Float64))
    shift_error < 1e-12 || throw(ArgumentError("LO shifts lose weighted adjointness: $shift_error"))
    weights = (energy_weights*radial_weights').*(spin_degeneracy/(2π))
    return LinkedProblem(collect(Float64,energy),collect(Float64,energy_weights),
        collect(Float64,radial_weights),weights,geo,convert(Array{ComplexF64,3},hamiltonian),
        convert(Matrix{Float64},hopping),convert(Array{ComplexF64,4},anchor),kernels,
        ne,nk,r,Int(block_width),nz,Float64(spin_degeneracy),Float64(current_scale),
        Float64(density_scale),Float64(target_dimless),shift_error)
end

function node_contract(ws,geometry)
    r,l,g = ntuple(j -> view(ws.local_green,:,:,j),3)
    radius = max(norm(r),floatmin(Float64))
    jump = norm((g-l)-(r-r'))/radius
    weighted = imag(tr(ws.mode_green[:,:,2]*geometry.center_metric))
    trace = imag(tr(l))
    number = abs(weighted-trace)/max(abs(trace),1e-300)
    # Sign of spectral component is assessed in the positive kinetic product;
    # checking the cancellation-prone difference of retarded blocks is separate.
    minimum_psd = zeros(3)
    for (i,p) in enumerate((im*(g-l),-im*l,im*g))
        ev = eigvals(Hermitian(p))
        minimum_psd[i] = min(0,first(ev)/max(maximum(abs,ev),1e-300))
    end
    return jump,number,minimum_psd
end

pack_fro_weights(rank) = vcat(ones(rank),fill(2.,rank^2-rank))

function weighted_packed_norm(matrix,weight,fro)
    return sum((fro*weight').*abs2.(matrix))
end

function source_residual(source,new,p,plan)
    d = pack_fro_weights(p.r)
    central = div(length(p.geo.offsets),2)*p.ne .+ (1:p.ne)
    h = view(plan.pv,central,:)
    output = Dict{String,Vector{Float64}}()
    temp = zeros(Float64,p.r^2)
    for (c,name) in enumerate(COLLISION_CHANNELS)
        diff,base = zeros(3),zeros(3)
        @views for k in 1:p.nk
            measure = p.weights[:,k]
            old_lg = source[:,:,:,k,c]
            new_lg = new[:,:,:,k,c]
            for component in 1:2
                diff[component+1] += weighted_packed_norm(new_lg[:,component,:]-old_lg[:,component,:],measure,d)
                base[component+1] += weighted_packed_norm(new_lg[:,component,:],measure,d)
            end
            gamma = new_lg[:,1,:]+new_lg[:,2,:]
            delta_gamma = gamma-(old_lg[:,1,:]+old_lg[:,2,:])
            pack_hermitian!(temp,p.const_by_channel[:,:,k,c])
            new_re = gamma*h'.+temp
            delta_re = delta_gamma*h'
            diff[1] += weighted_packed_norm(delta_re,measure,d)+weighted_packed_norm(delta_gamma,measure,d)/4
            base[1] += weighted_packed_norm(new_re,measure,d)+weighted_packed_norm(gamma,measure,d)/4
        end
        output[String(name)] = sqrt.(diff./max.(base,1e-300))
    end
    return output
end

function local_collision_balance(p,normalized,fresh)
    d = pack_fro_weights(p.r)
    output = Dict{String,Float64}()
    for (c,name) in enumerate(COLLISION_CHANNELS)
        diff,base = 0.,0.
        @inbounds for k in 1:p.nk, e in 1:p.ne
            gain,loss = 0.,0.
            for i in eachindex(d)
                gain += d[i]*fresh[i,1,e,k,c]*normalized[i,2,e,k]
                loss += d[i]*fresh[i,2,e,k,c]*normalized[i,1,e,k]
            end
            diff += p.weights[e,k]*(gain-loss)
            base += p.weights[e,k]*(abs(gain)+abs(loss))
        end
        output[String(name)] = abs(diff)/max(base,1e-300)
    end
    return output
end

function nodal_image(p,total,plan,anchor;previous=nothing,normalized=nothing,observe=_->nothing)
    (;r,ne,nk,nb,nz,geo) = p
    m = size(p.h,1)
    comp = p.density_scale
    j0 = p.current_scale
    stage_is_initial = isnothing(previous)
    packed = Array{Float64}(undef,r^2,2,ne,nk)
    ret = Array{ComplexF64}(undef,r,r,ne,nk)
    cavities = Array{ComplexF64}(undef,nb,nb,3,3,ne,nk)
    stored_response = stage_is_initial ? Array{ComplexF64}(undef,r,m,ne,nk) : nothing
    curr = zeros(nz+1,2)
    number,holes = 0.,0.
    metric = spectral_completeness(geo)
    base = centered_basis(geo)
    gr_defect = zeros(nk)
    max_spectral_identity,max_number,max_dyson,max_cavity = 0.,0.,0.,0.
    min_psd = zeros(3)
    ws = TailWorkspace(geo)
    rk_diff,rk_base = 0.,0.
    minus_i_lesser = zeros(ComplexF64,r,r)
    @views for k in 1:nk
        spectral = zeros(ComplexF64,r,r)
        source_ring = sample_packed_ring(plan,total[:,:,:,k],anchor[:,:,k])
        for e in 1:ne
            sample_node!(ws,geo,p.en[e],p.h[:,:,k],p.hops[:,k],source_ring[:,:,:,e,:])
            R,L,G = ws.local_green[:,:,1],ws.local_green[:,:,2],ws.local_green[:,:,3]
            ret[:,:,e,k] = R
            pack_hermitian!(packed[:,1,e,k],-im .* L)
            pack_hermitian!(packed[:,2,e,k], im .* G)
            parts = boundary_parts(ws,geo)
            cavities[:,:,:,:,e,k] = parts
            center = geo.center_modes
            edge = dropdims(sum(parts;dims=4);dims=4)
            s_r = inv(ws.d[center,center]-edge[:,:,1])
            s_l = s_r*(ws.mode_source[center,center,2]+edge[:,:,2])*s_r'
            s_g = s_r*(ws.mode_source[center,center,3]+edge[:,:,3])*s_r'
            max_cavity = max(max_cavity, norm(s_r-ws.mode_green[center,center,1])/max(norm(s_r),1e-300),
                norm(s_l-ws.mode_green[center,center,2])/max(norm(s_l),1e-300),
                norm(s_g-ws.mode_green[center,center,3])/max(norm(s_g),1e-300))
            max_dyson = max(max_dyson,norm(ws.d*ws.mode_green[:,:,1]-I)/sqrt(m))
            weight = p.weights[e,k]
            n,empty = imag(tr(L)),-imag(tr(G))
            number += weight*n
            holes += weight*empty
            curr .+= (j0*weight).*ws.bond_profile
            spectral .+= (im*p.we[e]/(2π)).*(R-R')
            jump,num,packed_psd = node_contract(ws,geo)
            max_spectral_identity = max(max_spectral_identity,jump)
            max_number = max(max_number,num)
            min_psd .= min.(min_psd,packed_psd)
            if stage_is_initial
                stored_response[:,:,e,k] .= base*ws.mode_green[:,:,1]
            else
                old_p = previous.response[:,:,e,k]
                candidate = -im .* (old_p*ws.mode_source[:,:,2]*old_p')
                unpack_hermitian!(minus_i_lesser,normalized[:,1,e,k])
                rk_diff += weight*sum(abs2,candidate-minus_i_lesser)
                rk_base += weight*sum(abs2,minus_i_lesser)
            end
        end
        gr_defect[k] = norm(spectral-metric)/norm(metric)
        observe((;k,spectral_defect=gr_defect[k]))
    end
    raw = Dict("raw_sheet_cm2"=>comp*number,"raw_holes_cm2"=>comp*holes,
        "raw_lambda"=>number/p.target_dimless,"raw_profile_A_cm2"=>curr[:,1],
        "raw_outward_left_right_A_cm2"=>[-curr[1,1],curr[end,1]],
        "spectral_weight_relative_by_ring"=>gr_defect,"kinetic_PSD_min_relative"=>min_psd,
        "Dyson_max_relative"=>max_dyson,"charge_metric_max_relative"=>max_number,
        "spectral_identity_max_relative"=>max_spectral_identity,"coupled_cavity_max_relative"=>max_cavity)
    if max(max_dyson,max_number,max_spectral_identity,max_cavity)>=1e-8 || minimum(min_psd)<-1e-9
        error("outer linked algebraic contracts failed: $raw")
    end
    return (;occupation=packed,retarded=ret,response=stored_response,cavities,raw,number,holes,curr,
        rK=stage_is_initial ? 0. : sqrt(rk_diff/rk_base))
end

function project_pair!(out, occupation, a, c)
    @views begin
        out[:,1,:,:] .= a.*occupation[:,1,:,:].+c.*occupation[:,2,:,:]
        out[:,2,:,:] .= (1-c).*occupation[:,2,:,:].+(1-a).*occupation[:,1,:,:]
    end
    return out
end

function prepare_occupation!(out,before,p)
    n,g,N = before.number,before.holes,p.target_dimless
    n>0 && g>0 || error("nonpositive raw occupation measure")
    a,c = min(1.,N/n),max(0.,(N-n)/g)
    c<=1 || error("No hole capacity for occupation projection")
    project_pair!(out,before.occupation,a,c)
    @views profile = a.*before.curr[:,1].-c.*before.curr[:,2]
    capture = sum(p.weights.*dropdims(sum(view(out,1:p.r,1,:,:);dims=1);dims=1))
    abs(capture/N-1)<1e-12 || error("normalization missed charge target")
    return Dict("physical_sheet_cm2"=>capture*p.density_scale,"a"=>a,"c"=>c,
        "outward_left_right_A_cm2"=>[-first(profile),last(profile)],"profile_A_cm2"=>profile)
end

function cavity_residuals(input,fresh,p)
    d,b = zeros(3,3),zeros(3,3)
    @views for k in 1:p.nk, e in 1:p.ne, side in 1:3, comp in 1:3
        w = p.weights[e,k]
        a = fresh[:,:,comp,side,e,k]
        d[comp,side] += w*sum(abs2,a-input[:,:,comp,side,e,k])
        b[comp,side] += w*sum(abs2,a)
    end
    return Dict(name=>sqrt.(d[:,i]./max.(b[:,i],1e-300)) for (i,name) in
        enumerate(("left","right","mixed")))
end


"""
    linked_fresh_map(problem, seed; projector=project_linked, selected_rows=false, ...)

Run one full unmixed map using rank² × LG × E × k × channel packed real sources.
Return the fresh packed sources, nodal Green functions and cavities, the adaptive
raw Green functions when enabled, and scalar/current-profile diagnostics.
`projector` is the ring projector used for both images. It receives the
`sampler` keyword to support the exact selected-row solve. Pass `nothing` to use
nodal observations alone. Set `projection_for_collision=false` to apply the
adaptive charge coefficients to nodal occupations before recomputing collisions.
The optional `observe` callback receives stage and ring or energy-tile progress.

The resonance projector defaults to order 12 on one subdivision. The
distance-tiered projector defaults to GL9/13/18 on two subdivisions.
"""
function linked_fresh_map(p::LinkedProblem, seed::AbstractArray{Float64,5};
        tile=8, projector=project_linked, projection_for_collision=true, selected_rows=false,
        order=projector === project_linked_economical ? 18 : 12,
        subdivisions=projector === project_linked_economical ? 2 : 1,
        grading=true, workers=Threads.nthreads(:default),
        lo_occupation_transfer=nothing, observe=_->nothing)
    size(seed) == (p.r^2,2,p.ne,p.nk,length(COLLISION_CHANNELS)) ||
        throw(DimensionMismatch("seed must be rank² × LG × E × k × channel"))
    tile > 0 && subdivisions > 0 && order >= 2 && workers > 0 ||
        throw(ArgumentError("invalid tile, subdivisions, quadrature order, or worker count"))
    projector === project_linked_economical && order != 18 &&
        throw(ArgumentError("distance-tiered linked projection uses order 18 near poles"))
    start = time()
    plan = CausalRows(p.en,shifted_queries(p.geo,p.en))
    anchor = dropdims(sum(p.const_by_channel;dims=4);dims=4)
    seed_total = dropdims(sum(seed;dims=5);dims=5)
    input = nodal_image(p,seed_total,plan,anchor;
        observe=record->observe((;stage=:input_nodal,record...)))
    projected_input = projected_image(p,seed_total,anchor,:input_projected;
        projector,selected_rows,order,subdivisions,grading,workers,observe)
    input_for_measure = isnothing(projected_input) ? input : projected_input
    normalized = similar(input.occupation)
    normalization = prepare_occupation!(normalized,input_for_measure,p)
    nodal_normalized = if isnothing(projected_input)
        normalized
    else
        project_pair!(similar(normalized),input.occupation,normalization["a"],normalization["c"])
    end
    collision_occupation = projection_for_collision ? normalized : nodal_normalized
    source = Array{Float64}(undef,size(seed))
    scatter_packed!(source,collision_occupation,p.v;tile,lo_occupation_transfer,
        observe=record->observe((;stage=:collision,record...)))
    balances = local_collision_balance(p,collision_occupation,source)
    fresh_total = dropdims(sum(source;dims=5);dims=5)
    fresh = nodal_image(p,fresh_total,plan,anchor;previous=input,normalized=nodal_normalized,
        observe=record->observe((;stage=:fresh_nodal,record...)))
    projected_fresh = projected_image(p,fresh_total,anchor,:fresh_projected;
        projector,selected_rows,order,subdivisions,grading,workers,observe)
    residuals = source_residual(seed,source,p,plan)
    boundary = cavity_residuals(input.cavities,fresh.cavities,p)
    summary = Dict{String,Any}(
        "voltage_meV_per_period"=>p.geo.voltage*100,"dimensions"=>[p.ne,p.nk,p.r,p.nb],
        "physical_E_measure_meV"=>100sum(p.we),"physical_radial_measure"=>sum(p.wk),
        "LO_weighted_adjoint_relative"=>p.shift_error,
        "LO_occupation_transfer"=>isnothing(lo_occupation_transfer) ?
            "standard Bose-weighted adjoint pair" : "explicit custom LO occupation matrices",
        "input"=>input.raw,"input_pair_projection"=>normalization,
        "fresh"=>fresh.raw,"weighted_rK"=>fresh.rK,"weighted_rSigma_channels_RLG"=>residuals,
        "weighted_rSigma_boundaries_RLG"=>boundary,"fresh_channel_balance"=>balances,
        "elapsed_seconds_including_setup"=>time()-start)
    if !isnothing(projector)
        summary["projected_input"] = projected_input.raw
        summary["projected_fresh"] = projected_fresh.raw
        summary["rK_measure_convention"] =
            "nodal fresh candidate on input retarded, with adaptive projection coefficients"
        summary["collision_occupation"] = projection_for_collision ?
            "positive adaptive hat projection" : "nodal occupation with adaptive charge coefficients"
    end
    fresh_state = (;occupation=fresh.occupation,retarded=fresh.retarded,cavities=fresh.cavities)
    adaptive_state = isnothing(projected_fresh) ? nothing :
        (;occupation=projected_fresh.occupation,retarded=projected_fresh.retarded,
            seam_lesser=projected_fresh.seam_lesser)
    return (;source,fresh=fresh_state,projected=adaptive_state,summary)
end

function projected_image(p,total,anchor,stage;projector,selected_rows,order,subdivisions,
        grading,workers,observe)
    isnothing(projector) && return nothing
    occupation = similar(total)
    retarded = Array{ComplexF64}(undef,p.r,p.r,p.ne,p.nk)
    return linked_projected_stage!(occupation,p,total,anchor;retarded,projector,selected_rows,
        order,subdivisions,grading,workers,observe=record->observe((;stage,record...)))
end
