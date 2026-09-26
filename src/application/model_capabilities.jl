"""
    model_capabilities(problem, algorithms)

Machine-readable model scope. Numerical convergence certifies a fixed point
of the declared discrete equations, independently of physical model adequacy.
These notices describe model assumptions; they are not failed convergence gates.
"""
function model_capabilities(problem::NEGFProblem, algorithms::AlgorithmOptions)
    limits = Dict{String,Any}[]
    note(code, message) = push!(limits, Dict("code"=>code, "message"=>message))
    localization = problem.basis.localization
    if localization in (:pzp,)
        note(
            "BASIS_CELL_CLIPPING",
            "Multi-period states are clipped and reorthogonalized in one cell; inspect independent full-window and BDD dispersion references.",
        )
    elseif localization === :none
        note(
            "TRUNCATED_CELL_EIGENSTATES",
            "A finite cell eigenstate subspace is used; omitted states require an independent basis convergence study.",
        )
    end
    note(
        "FINITE_SPATIAL_GRID",
        "The real-space BDD grid and energy/momentum windows require separate refinement studies.",
    )
    if algorithms.embedding === :full_cell_resolvent
        note(
            "BULK_GREEN_FEEDBACK",
            "The legacy neighbouring bulk Green feedback is not an exact surface or cavity self-energy.",
        )
    else
        note(
            "FINITE_OPEN_CHAIN",
            "Cavity elimination is exact for the finite local-self-energy chain; increase neighbours per side and energy window to test truncation.",
        )
    end
    note(
        "LOCAL_SCATTERING_VERTICES",
        "Scattering self-energies are local to a cell; full-window reference bases are not a completed nonlocal transport implementation.",
    )
    if problem.scattering.LO
        note(
            "LO_LEGACY_REGULATOR",
            "LO uses the existing 1/(Q^2+qs^2) regulator, not the distinct statically screened Frohlich factor Q^2/(Q^2+qs^2)^2.",
        )
        note(
            "LO_QUADRATURE_UNCERTIFIED",
            "LO angular and longitudinal quadratures require their own convergence study.",
        )
    end
    if problem.scattering.impurity && algorithms.kernel_build !== :adaptive_direct
        note(
            "IMPURITY_ANGLE_FINITE",
            "Finite-angle agreement or zero interpolation error does not certify the continuum forward-scattering integral.",
        )
    end
    if algorithms.hilbert in (:direct, :fft)
        note(
            "LEGACY_HILBERT_DISCRETIZATION",
            "The omitted-diagonal Hilbert sum is retained as a comparison discretization; compare product integration and enlarge the energy window.",
        )
    end
    if problem.models.electron_electron.mode === :none
        note("NO_ELECTRON_ELECTRON", "Electron-electron scattering is not included.")
    else
        note(
            "PRESCRIBED_SPPA_BATH",
            "Single-q SPPA exchanges energy with a prescribed plasmon bath; it is not a self-consistent conserving GW electron interaction.",
        )
    end
    if algorithms.occupation_normalization === :scalar_lesser
        note(
            "LEGACY_LESSER_NORMALIZATION",
            "Unilateral lesser normalization is retained only as a controlled legacy comparison; it can violate empty-state positivity before convergence.",
        )
    end
    note(
        "BARE_BUBBLE_OPTICS",
        "Optical response has no conserving vertex correction and is not a certified absolute experimental gain prediction.",
    )
    return Dict{String,Any}(
        "schema"=>"qcl-negf-model-capabilities-v1",
        "localization"=>String(localization),
        "embedding"=>String(algorithms.embedding),
        "embedding_neighbours_per_side"=>algorithms.embedding_periods,
        "energy_shift"=>String(algorithms.energy_shift),
        "hilbert"=>String(algorithms.hilbert),
        "kernel_build"=>String(algorithms.kernel_build),
        "occupation_normalization"=>String(algorithms.occupation_normalization),
        "quantitative_physics_certified"=>false,
        "reference_only_bases"=>["pzp_tails", "bloch_wannier", "wannier_stark"],
        "known_limitations"=>limits,
    )
end
