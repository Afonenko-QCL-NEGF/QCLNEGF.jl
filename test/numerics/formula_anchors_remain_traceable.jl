module Suite_T103
include("../support/common.jl")
include("../support/traceability.jl")

@testset "Formula anchors remain traceable" begin
    root = normpath(joinpath(TEST_ROOT, ".."))
    docs = join(
        read(joinpath(root, "docs", "src", "theory", file), String) for
        file in readdir(joinpath(root, "docs", "src", "theory"))
    )
    source_files = String[]
    for (directory, _, files) in walkdir(joinpath(root, "src"))
        append!(
            source_files,
            joinpath.(directory, filter(file -> endswith(file, ".jl"), files)),
        )
    end
    sources = join(read(file, String) for file in sort!(source_files))
    contracts = [
        ("eq-bdd-matrix", "build_bdd_hamiltonian"),
        ("eq-pzp-lowdin", "build_localized_basis"),
        ("eq-reduced-hamiltonian", "project_hamiltonians"),
        ("eq-energy-shift", "build_shift_matrix"),
        ("eq-embedding", "embedding_self_energy"),
        ("eq-dyson-retarded", "retarded_green"),
        ("eq-keldysh", "keldysh_green"),
        ("eq-number-normalization", "number_functional"),
        ("eq-realspace-density", "electron_density"),
        ("eq-kernel-rescale", "build_kernels"),
        ("eq-impurity-kernel", "impurity_kernel"),
        ("eq-ifr-kernel", "interface_roughness_kernel"),
        ("eq-acoustic-kernel", "acoustic_kernel"),
        ("eq-alloy-kernel", "alloy_kernel"),
        ("eq-lo-kernel", "lo_phonon_kernel"),
        ("eq-static-contraction", "_static_contraction"),
        ("eq-lo-contraction", "_lo_contraction"),
        ("eq-direct-hilbert", "direct_hilbert_transform"),
        ("eq-retarded-selfenergy", "retarded_self_energy"),
        ("eq-scba-map", "solve_scba"),
        ("eq-poisson-matrix", "build_poisson_matrix"),
        ("eq-bordered-poisson", "solve_periodic_poisson"),
        ("eq-outer-map", "solve"),
        ("eq-sheet-density", "sheet_density_matrix"),
        ("eq-boundary-flux", "boundary_flux"),
        ("eq-resolved-current", "energy_resolved_current"),
        ("eq-collision-balance", "collision_balance"),
        ("eq-power-balance", "power_balance"),
        ("eq-effective-levels", "effective_levels"),
        ("eq-validation-residuals", "validate_solution"),
    ]
    for (anchor, binding) in contracts
        @test occursin("(@id $anchor)", docs)
        @test occursin(binding, sources)
    end
end

end # independent suite
