function observation(geo::TailGeometry)
    lo, hi = first(geo.monitor), last(geo.monitor)
    return vcat(centered_basis(geo),geo.modes[[lo,lo+1,hi-1,hi],:])
end

"""Use D′ with the central and both seam coordinate observations as right-hand sides."""
function observed_sampler(geo::TailGeometry; audit=nothing)
    obs = observation(geo)
    nobs, n = size(obs)
    rank = nobs-4
    right_hand_side = copy(obs')
    solved = similar(right_hand_side)
    observed = zeros(ComplexF64,nobs,nobs)
    scratch = zeros(ComplexF64,nobs,n)

    return function sample_observed!(ws::TailWorkspace,geo,energies,h,hopping,source,rows)
        nenergy = length(energies)
        fields = (source.Σᴿ,source.Σˡ,source.Σᵍ)
        result = ntuple(_->Array{ComplexF64}(undef,nenergy,rank,rank),3)
        seams = zeros(Float64,3,2,nenergy)
        seam_lesser = zeros(2,nenergy)
        length(rows)==length(geo.offsets) || throw(DimensionMismatch("shifted periods don't match geometry"))
        @views for e in eachindex(energies)
            fill!(ws.triplets,0)
            for (shift,plan) in enumerate(rows), (row,weight) in plan[e], comp in 1:3
                ws.triplets[:,:,comp,shift] .+= weight .* fields[comp][row,1,:,:]
            end
            assemble_sources!(ws,geo,ws.triplets)
            ws.d .= energies[e].*geo.metric .-h.-ws.mode_source[:,:,1]
            copyto!(ws.factor,ws.d')
            factor = lu!(ws.factor)
            copyto!(solved,right_hand_side)
            ldiv!(factor,solved)
            mul!(observed,solved',right_hand_side)
            result[1][e,:,:] .= observed[1:rank,1:rank]
            for comp in 2:3
                mul!(scratch,solved',ws.mode_source[:,:,comp])
                mul!(observed,scratch,solved)
                result[comp][e,:,:] .= observed[1:rank,1:rank]
                left = 2first(hopping)*real(observed[rank+2,rank+1])
                right = 2last(hopping)*real(observed[rank+4,rank+3])
                if comp==2
                    seam_lesser[1,e] = real(observed[rank+2,rank+1])
                    seam_lesser[2,e] = real(observed[rank+4,rank+3])
                    seams[1:2,1,e] .= right
                    seams[1:2,2,e] .= -left
                else
                    seams[3,1,e] = -right
                    seams[3,2,e] = left
                end
            end
            isnothing(audit) || audit(e,ws.d,solved,right_hand_side)
        end
        return (;gr=result[1],gl=result[2],gg=result[3],seams,seam_lesser)
    end
end
