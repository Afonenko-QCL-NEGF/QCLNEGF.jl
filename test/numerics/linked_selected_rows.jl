using LinearAlgebra
using Test
using QCLNEGF.LinkedTransport: LinkedProblem, LocalCollisionKernels, linked_fresh_map

@testset "Full and selected linked solves agree on one weighted map" begin
    energy = [-4.0, -3.4, -2.0, -1.0, 0.0, 0.7, 2.0, 3.2, 4.0]
    widths = diff(energy)
    weights_e = [widths[1]; widths[1:end-1] + widths[2:end]; widths[end]] ./ 2
    weights_k = [0.25, 0.75]
    frame = reshape(ComplexF64.(normalize([1.0, 0.6])), 2, 1)
    modes = zeros(ComplexF64, 6, 3)
    for period in 1:3
        modes[2period-1:2period, period] .= vec(frame)
    end
    nk, ne = length(weights_k), length(energy)
    h = cat(
        Matrix(SymTridiagonal([0.2, 0.3, 0.4], [0.08, 0.12])),
        Matrix(SymTridiagonal([0.3, 0.4, 0.6], [0.10, 0.14]));
        dims=3,
    )
    lo_transfer = Matrix{Float64}(I, ne, ne)
    kernels = LocalCollisionKernels(
        frame, weights_k, [0.05, 0.01], [0.005, 0.008], [0.6 0.4],
        fill(0.012, nk, nk), fill(0.018, 2, 2, nk, nk), fill(0.009, 2, 2, nk, nk),
        lo_transfer, lo_transfer, 0.12,
    )
    hops = [0.12 0.16; 0.10 0.11; 0.14 0.18]
    anchor = zeros(ComplexF64, 1, 1, nk, 5)
    scale = (; spin_degeneracy=2, current_scale=1.0, density_scale=1.0, target_dimless=0.1)
    p = LinkedProblem(energy, weights_e, weights_k, modes, 1, [-1, 0, 1], 0.6, 1,
        h, hops, anchor, kernels; scale...)
    seed = zeros(1, 2, ne, nk, 5)
    seed[1, 1, :, :, 1] .= [0.085 0.12]
    seed[1, 2, :, :, 1] .= [0.25 0.16]

    full = linked_fresh_map(p, seed; order=12, subdivisions=2, workers=1)
    selected = linked_fresh_map(p, seed; order=12, subdivisions=2, workers=1, selected_rows=true)
    @test norm(full.summary["projected_fresh"]["raw_profile_A_cm2"]) > 1e-5
    @test all(norm(selectdim(full.source, 5, channel)) > 1e-5 for channel in 1:5)
    for (a, b) in (
        (full.source, selected.source),
        (full.projected.occupation, selected.projected.occupation),
        (full.projected.retarded, selected.projected.retarded),
        (full.projected.seam_lesser, selected.projected.seam_lesser),
        (full.summary["weighted_rK"], selected.summary["weighted_rK"]),
    )
        # Roundoff allowance on a small Float64 solve in native dimensionless units.
        @test isapprox(a, b; rtol=2e-10, atol=1e-12)
    end
    for stage in ("projected_input", "projected_fresh")
        for field in ("raw_sheet_cm2", "raw_lambda", "raw_profile_A_cm2",
                      "raw_greater_profile_A_cm2")
            @test isapprox(full.summary[stage][field], selected.summary[stage][field];
                rtol=2e-10, atol=1e-12)
        end
    end
    @test_throws ArgumentError LinkedProblem(energy, weights_e, reverse(weights_k), modes,
        1, [-1, 0, 1], 0.6, 1, h, hops, anchor, kernels; scale...)
end
