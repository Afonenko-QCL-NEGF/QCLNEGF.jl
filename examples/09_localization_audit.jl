#!/usr/bin/env julia
# Small independent basis study. No engine, web, HDF5 or service account.
# julia --project=. examples/09_localization_audit.jl /tmp/reference2019-basis-audit
using QCLNEGF
using LinearAlgebra
using Unitful

output = abspath(isempty(ARGS) ? "results/localization-audit" : only(ARGS))
mkpath(output)
p, s = reference_parameters(), ScaleSystem()

open(joinpath(output, "localization_summary.csv"), "w") do table
    println(
        table,
        "method,spatial_points_per_period,states_per_period,periods,orthogonality_residual,subspace_residual,maximum_spectral_error_eV,maximum_central_tail_weight,minimum_gauge_overlap",
    )
    for Nz in (24, 48), Nb in (3, 5), P in (1, 2)
        n = NumericalParameters(
            N_z = Nz,
            N_b = Nb,
            P_basis = P,
            E_min = -0.2u"eV",
            E_max = 0.5u"eV",
            N_E = 31,
            M_E = 0.08u"eV",
            k_max = 0.4u"nm^-1",
            N_k = 3,
            N_φ = 8,
            qz_max = 5u"nm^-1",
            N_qz = 9,
            η_seed = 2u"meV",
        )
        grids = build_grids(p, n, s)
        profiles = build_profiles(p, n, s, grids)
        for method in (:pzp_tails, :bloch_wannier, :wannier_stark)
            reference = QCLNEGF.build_multiperiod_reference_basis(
                p,
                n,
                s,
                grids,
                profiles;
                localization = method,
            )
            spectral_error =
                maximum(
                    abs,
                    eigvals(Hermitian(reference.hamiltonian)) -
                    reference.reference_energies,
                ) * s.E₀_eV
            println(
                table,
                join(
                    (
                        method,
                        Nz,
                        Nb,
                        reference.periods,
                        norm(reference.overlap-I),
                        reference.subspace_residual,
                        spectral_error,
                        maximum(reference.central_tail_weights),
                        something(reference.gauge_minimum_overlap, ""),
                    ),
                    ',',
                ),
            )
            prefix = "$(method)-Nz$(Nz)-Nb$(Nb)-P$(reference.periods)"
            open(joinpath(output, prefix * "-envelopes.csv"), "w") do envelopes
                println(envelopes, "z_nm,state,real_envelope,imag_envelope")
                for state in reference.central_columns,
                    g in eachindex(reference.coordinates)

                    wave = reference.envelopes[g, state]
                    println(
                        envelopes,
                        join(
                            (
                                reference.coordinates[g]*s.L₀_m*1e9,
                                state,
                                real(wave),
                                imag(wave),
                            ),
                            ',',
                        ),
                    )
                end
            end
        end
        for method in (:legacy_pzp, :none)
            basis = build_basis(p, n, s, grids, profiles; localization = method)
            check = QCLNEGF.basis_dispersion_diagnostic(basis, profiles, grids, s)
            # A dispersion error compares with full BDD; it is not equivalent
            # to the invariant-subspace residual used in the full-window rows.
            println(
                table,
                join(
                    (
                        method,
                        Nz,
                        Nb,
                        2P+1,
                        basis.r_orth,
                        "",
                        check.maximum_error_eV,
                        "",
                        "",
                    ),
                    ',',
                ),
            )
        end
        flush(table)
    end
end
println("Basis study written to ", output)
println("This is an operator/subspace diagnostic, not a converged transport study.")
