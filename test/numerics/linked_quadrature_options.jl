using LinearAlgebra
using Test
import QCLNEGF
using QCLNEGF.LinkedTransport: LinkedProblem, LocalCollisionKernels, linked_fresh_map, project_linked_economical

@testset "Linked distance projector uses its own defaults" begin
    energy = collect(-4.0:4.0)
    weights_e = [0.5; ones(7); 0.5]
    weights_k = [1.0]
    coordinate_frame = reshape(ComplexF64[1,0],2,1)
    modes = zeros(ComplexF64,6,3)
    for (period, row) in enumerate(1:2:5)
        modes[row,period] = 1
    end
    nk, ne = length(weights_k), length(energy)
    h = reshape(Matrix{ComplexF64}(Diagonal([0.2,0.3,0.4])),3,3,nk)
    lo_transfer = Matrix{Float64}(I,ne,ne)
    kernels = LocalCollisionKernels(coordinate_frame,weights_k,[0.05,0.0],zeros(2),
        zeros(1,2),zeros(nk,nk),zeros(ComplexF64,2,2,nk,nk),zeros(2,2,nk,nk),
        lo_transfer,lo_transfer,0.0)
    p = LinkedProblem(energy,weights_e,weights_k,modes,1,[-1,0,1],1.0,1,
        h,zeros(3,nk),zeros(ComplexF64,1,1,nk,5),kernels;
        spin_degeneracy=2,current_scale=1,density_scale=1,target_dimless=0.1)
    seed = zeros(1,2,ne,nk,5)
    seed[1,1,:,:,1] .= 0.1
    seed[1,2,:,:,1] .= 0.2

    projector = project_linked_economical
    implicit = linked_fresh_map(p,seed;projector,workers=1)
    explicit = linked_fresh_map(p,seed;projector,order=18,subdivisions=2,workers=1)
    @test implicit.source == explicit.source
    @test implicit.projected.occupation == explicit.projected.occupation
    @test implicit.summary["projected_fresh"]["sample_count_by_ring"] ==
        explicit.summary["projected_fresh"]["sample_count_by_ring"]
    started = Ref(false)
    @test_throws ArgumentError linked_fresh_map(p,seed;projector,order=12,
        observe=_->(started[]=true))
    @test !started[]
end
