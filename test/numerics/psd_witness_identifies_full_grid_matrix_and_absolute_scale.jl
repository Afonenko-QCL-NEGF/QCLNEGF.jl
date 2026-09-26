module Suite_T097
include("../support/common.jl")

@testset "PSD witness identifies full-grid matrix and absolute scale" begin
    NE, Nk, Nb=3, 2, 2
    ε=collect(range(-1.0, 1.0; length = NE))
    h=zeros(ComplexF64, Nk, Nb, Nb)
    GR=zeros(ComplexF64, NE, Nk, Nb, Nb)
    GL=similar(GR)
    GG=similar(GR)
    A=similar(GR)
    total=SelfEnergyFamily(
        zeros(ComplexF64, size(GR)),
        zeros(ComplexF64, size(GR)),
        zeros(ComplexF64, size(GR)),
    )
    for e = 1:NE, m = 1:Nk
        GR[e, m, :, :] .= -im .* Matrix{ComplexF64}(I, Nb, Nb)
        GL[e, m, :, :] .= im .* Matrix{ComplexF64}(I, Nb, Nb)
        GG[e, m, :, :] .= -im .* Matrix{ComplexF64}(I, Nb, Nb)
        A[e, m, :, :] .= 2 .* Matrix{ComplexF64}(I, Nb, Nb)
        total.Σᴿ[e, m, :, :] .= -im .* Matrix{ComplexF64}(I, Nb, Nb)
        total.Σˡ[e, m, :, :] .= im .* Matrix{ComplexF64}(I, Nb, Nb)
        total.Σᵍ[e, m, :, :] .= -im .* Matrix{ComplexF64}(I, Nb, Nb)
    end
    green=GreenState(GR, GL, GG, A, ones(NE, Nk), ones(NE, Nk))
    # Unoccupied spectral density iG> has a negative eigenvalue at E3,k2.
    GG[3, 2, 1, 1]=2im
    original=copy(GG)
    reference=QCLNEGF._green_positivity_witnesses(green, total)
    witness=reference[argmax([item.ratio for item in reference])]
    @test witness.matrix_kind===:unoccupied
    @test (witness.energy_index, witness.momentum_index)==(3, 2)
    @test witness.minimum_eigenvalue==-2.0
    @test witness.block_norm==2.0
    @test witness.absolute_defect == 2.0
    @test witness.relative_defect == 1.0
    @test witness.backward_error > 0
    @test witness.ratio ≈ (2.0 - witness.backward_error) / 2.0
    @test witness.ratio > 1e-10
    @test GG==original
    for chunk_size in (1, 4, 6)
        workspace=ProductionResidualWorkspace(NE, Nk, Nb; chunk_size)
        result=production_residual_suite(ε, h, green, total, total; workspace)
        production=QCLNEGF._production_positivity_witness(workspace)
        @test production==witness
        @test production.ratio==result.r_PSD
        # Workspace reuse resets the witness rather than leaking an old maximum.
        GG[3, 2, 1, 1]=-im
        clean=production_residual_suite(ε, h, green, total, total; workspace)
        @test clean.r_PSD==0.0
        @test QCLNEGF._production_positivity_witness(workspace).ratio==0.0
        GG[3, 2, 1, 1]=2im
    end
    # Dimensionless norm carries the missing physical magnitude: a ratio near
    # one can arise from a tiny negative block as well as a macroscopic one.
    tiny=QCLNEGF._positivity_witness(:broadening, 1, 1, [-1e-7, -2e-8])
    @test tiny.ratio>0.99999
    @test tiny.block_norm==1e-7
    @test tiny.minimum_eigenvalue==-1e-7
    bad=copy(GG)
    bad[3, 2, 1, 1]=NaN+0im
    invalid=GreenState(GR, GL, bad, A, ones(NE, Nk), ones(NE, Nk))
    failed=QCLNEGF._green_positivity_witnesses(invalid, total)[3]
    @test failed.matrix_kind===:unoccupied
    @test (failed.energy_index, failed.momentum_index)==(3, 2)
    @test failed.ratio==Inf
end

end # independent suite
