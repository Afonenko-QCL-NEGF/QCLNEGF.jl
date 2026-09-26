module Suite_T088
include("../support/common.jl")
include("../support/production_residuals.jl")

@testset "Production self-energy residual reuses bounded storage" begin
    _, _, _, old, candidate =
        _random_residual_fixture(NE = 17, Nk = 4, Nb = 3, seed = 0x51CBA)
    workspace = ProductionResidualWorkspace(17, 4, 3; chunk_size = 5)
    result = production_selfenergy_residual(old, candidate; workspace)
    reference = QCLNEGF._selfenergy_residual(old, candidate)
    @test result ≈ reference rtol=2e-14 atol=2e-14
    @test QCLNEGF._selfenergy_residual_production(old, candidate; workspace) == result
    @test production_selfenergy_residual(candidate, candidate; workspace) == 0.0
    @test production_selfenergy_residual(old, candidate; workspace) == result

    # Warm the compiled path before measuring.  The bound excludes any
    # four-dimensional `candidate-old` temporary and remains independent of
    # the number of matrix entries; task bookkeeping may allocate a few KiB.
    production_selfenergy_residual(old, candidate; workspace)
    bytes = @allocated production_selfenergy_residual(old, candidate; workspace)
    @test bytes < 32768 + 4096Threads.nthreads()
    measured = Int[]
    for sizeE in (17, 65, 257)
        _, _, _, previous, next = _random_residual_fixture(NE = sizeE, Nk = 4, Nb = 3)
        local_workspace = ProductionResidualWorkspace(sizeE, 4, 3; chunk_size = 5)
        production_selfenergy_residual(previous, next; workspace = local_workspace)
        push!(
            measured,
            @allocated production_selfenergy_residual(
                previous,
                next;
                workspace = local_workspace,
            )
        )
    end
    @test maximum(measured) - minimum(measured) <= 8192
end

end # independent suite
