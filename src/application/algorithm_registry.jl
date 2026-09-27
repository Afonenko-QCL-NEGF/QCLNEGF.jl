"""
One auditable row in [`optimization_catalog`](@ref).  `impact` is one of
`:physics_preserving`, `:controlled_numerical`, `:physical_model`, or
`:research_only`; `status` is `:implemented` or `:documented`.

See [Optimization decision tree](@ref optimization-decision-tree).
"""
struct OptimizationDescriptor
    id::Symbol
    label::String
    impact::Symbol
    status::Symbol
    summary::String
    reference::String
end

const _OPTIMIZATION_CATALOG = OptimizationDescriptor[
    OptimizationDescriptor(
        :product_integration_hilbert,
        "Piecewise-linear PV integral",
        :controlled_numerical,
        :implemented,
        "Analytic principal-value product integration with explicit finite-window taper.",
        "QCLNEGF independent causal-function tests",
    ),
    OptimizationDescriptor(
        :adaptive_impurity_quadrature,
        "Adaptive impurity angular quadrature",
        :controlled_numerical,
        :implemented,
        "Refines the narrow forward-scattering peak.",
        "QCLNEGF angular-refinement tests",
    ),
    OptimizationDescriptor(
        :conservative_shift_pair,
        "Quadrature-adjoint energy shift pair",
        :controlled_numerical,
        :implemented,
        "Opposite shifts obey the weighted discrete adjoint contract.",
        "QCLNEGF detailed-balance operator tests",
    ),
    OptimizationDescriptor(
        :finite_chain_embedding,
        "Finite chain cavity embedding",
        :physical_model,
        :implemented,
        "Retarded/lesser/greater cavity embedding with explicit period truncation.",
        "QCLNEGF finite-chain reference tests",
    ),
    OptimizationDescriptor(
        :sparse_energy_shift,
        "Sparse interpolation plan",
        :physics_preserving,
        :implemented,
        "O(N_E) two-point energy shifts replacing the dense interpolation matrix.",
        "QCLNEGF discrete field-periodic equations",
    ),
    OptimizationDescriptor(
        :fft_hilbert,
        "Zero-padded FFT Hilbert transform",
        :physics_preserving,
        :implemented,
        "Linear convolution for the same discrete principal-value sum; only summation order changes.",
        "QCLNEGF direct/FFT roundoff contract",
    ),
    OptimizationDescriptor(
        :dense_blas,
        "Compound-index BLAS contraction",
        :physics_preserving,
        :implemented,
        "Reorders the literal six-index contraction into chunked matrix products.",
        "QCLNEGF production backend",
    ),
    OptimizationDescriptor(
        :threaded_blocks,
        "Threaded independent blocks",
        :physics_preserving,
        :implemented,
        "Static disjoint partitions of Dyson, Keldysh, embedding, residual, and contraction work.",
        "Julia task/thread execution; unchanged discrete equations",
    ),
    OptimizationDescriptor(
        :momentum_independent_factorization,
        "Exact momentum-independent factorization",
        :physics_preserving,
        :implemented,
        "Integrates k-prime before the basis contraction when the tensor is exactly independent of both momenta.",
        "Algebraic factorization of the six-index kernel",
    ),
    OptimizationDescriptor(
        :anderson_mixing,
        "Anderson acceleration",
        :physics_preserving,
        :implemented,
        "Accelerates the same fixed point; convergence still requires the unaccelerated SCBA residuals.",
        "D. G. Anderson, J. ACM 12, 547 (1965)",
    ),
    OptimizationDescriptor(
        :tabulated_q_kernel,
        "Adaptive q lookup",
        :controlled_numerical,
        :implemented,
        "Interpolates microscopic q-dependent blocks with measured off-grid residuals.",
        "QCLNEGF production kernel diagnostics",
    ),
    OptimizationDescriptor(
        :low_rank_kernel,
        "Truncated-SVD kernel operator",
        :controlled_numerical,
        :implemented,
        "Compresses a flattened scattering map with a stored Frobenius-tail certificate.",
        "L. Zeng et al., arXiv:1304.0316; implementation is operator SVD, not global LRA",
    ),
    OptimizationDescriptor(
        :drop_retarded_real_part,
        "Omit Kramers-Kronig term in Sigma^R",
        :physical_model,
        :implemented,
        "Omits the dispersive Kramers-Kronig contribution, retains -i Gamma/2, and can move resonances and gain peaks.",
        "Lemus, Charles, and Kubis, arXiv:2003.09536",
    ),
    OptimizationDescriptor(
        :diagonal_self_energy,
        "Diagonal self-energy",
        :physical_model,
        :implemented,
        "Zeros off-diagonal basis coherences after each scattering candidate.",
        "Kubis and Vogl, Phys. Rev. B 83, 195304 (2011)",
    ),
    OptimizationDescriptor(
        :momentum_averaging,
        "Radial-momentum averaged self-energy",
        :physical_model,
        :implemented,
        "Replaces each k-resolved candidate by its normalized radial quadrature average.",
        "Kubis and Vogl, Phys. Rev. B 83, 195304 (2011)",
    ),
    OptimizationDescriptor(
        :recursive_green_function,
        "Recursive/selected Green function",
        :research_only,
        :documented,
        "Useful for long block-tridiagonal real-space devices; the current localized reference design basis is dense and small.",
        "Svizhenko et al., J. Appl. Phys. 91, 2343 (2002)",
    ),
    OptimizationDescriptor(
        :selected_inversion,
        "Selected inversion",
        :research_only,
        :documented,
        "Computes inverse entries associated with sparse-factor fill without materializing the full inverse; it needs a large sparse Dyson operator.",
        "Lin et al., ACM TOMS 37, 40 (2011)",
    ),
    OptimizationDescriptor(
        :find,
        "FIND nested dissection",
        :research_only,
        :documented,
        "Uses separator trees to obtain selected Green-function entries in wide sparse devices; the current small dense state basis has no such tree.",
        "Li et al., J. Comput. Phys. 227, 9408 (2008)",
    ),
    OptimizationDescriptor(
        :contact_block_reduction,
        "Contact block reduction",
        :research_only,
        :documented,
        "Exploits low-rank open-contact self-energies in ballistic devices and does not directly apply to the closed field-periodic SCBA cell.",
        "D. Mamaluy, D. Vasileska, M. Sabathil, T. Zibold, and P. Vogl, Phys. Rev. B 71, 245321 (2005)",
    ),
    OptimizationDescriptor(
        :mode_space,
        "Localized mode/Wannier space",
        :research_only,
        :documented,
        "Can reduce Dyson rank but requires projected four-index interaction tensors and new validation.",
        "Lee and Wacker, arXiv:cond-mat/0109163; Lemus et al., arXiv:2003.09536",
    ),
    OptimizationDescriptor(
        :global_lra,
        "Global low-rank approximation",
        :research_only,
        :documented,
        "Projects the whole NEGF problem and requires controlled back-transformation of local observables.",
        "L. Zeng et al., arXiv:1304.0316",
    ),
    OptimizationDescriptor(
        :model_order_reduction,
        "Projection/model-order reduction",
        :research_only,
        :documented,
        "Builds a reduced energy/parameter response space; causality, conservation, and local-observable back-transforms require new contracts.",
        "Huang et al., IEEE TED 60, 2111 (2013); Chen et al., JCP 286, 49 (2015)",
    ),
    OptimizationDescriptor(
        :broyden_mixing,
        "Broyden quasi-Newton mixing",
        :research_only,
        :documented,
        "Updates an approximate inverse Jacobian of the nonlinear fixed-point residual; Anderson is the implemented accelerated mixer.",
        "C. G. Broyden, Math. Comp. 19, 577 (1965)",
    ),
    OptimizationDescriptor(
        :buttiker_probes,
        "Multi-scattering Buettiker model",
        :research_only,
        :documented,
        "A faster quasi-equilibrium lesser-self-energy model, not the reference SCBA equations.",
        "Greck et al., Opt. Express 23, 6587 (2015)",
    ),
    OptimizationDescriptor(
        :density_matrix,
        "Density-matrix surrogate",
        :research_only,
        :documented,
        "A distinct reduced transport model requiring its own validity domain.",
        "Soleimanikahnoj et al., arXiv:1710.08870",
    ),
    OptimizationDescriptor(
        :stochastic_rnegf,
        "Stochastic rNEGF",
        :research_only,
        :documented,
        "Trades deterministic work for sampling error; expected self-averaging is weak in this effective 1D model.",
        "Zhang et al., Phys. Rev. B 110, 155430 (2024)",
    ),
    OptimizationDescriptor(
        :mixed_precision,
        "Mixed precision with iterative refinement",
        :research_only,
        :documented,
        "Can preserve a high-precision residual only under conditioning and fallback gates; unconditional FP32 is not equivalent.",
        "Carson and Higham, SIAM J. Sci. Comput. 40, A817 (2018)",
    ),
    OptimizationDescriptor(
        :gpu_batched_backend,
        "Batched GPU backend",
        :research_only,
        :documented,
        "Would batch many energy-momentum contractions and solves while keeping data resident; isolated tiny Dyson blocks are insufficient.",
        "Sawant et al., npj Comput. Mater. 11, 110 (2025)",
    ),
    OptimizationDescriptor(
        :distributed_energy_momentum,
        "Distributed energy/momentum decomposition",
        :research_only,
        :documented,
        "Requires explicit shifted-energy halos, momentum reductions, and reproducible collectives for SCBA coupling.",
        "Steiger et al., IEEE Trans. Nanotechnol. 10, 1464 (2011)",
    ),
]

"""
    optimization_catalog()

Return a copy of the optimization taxonomy used by the documentation and
machine-readable run manifest.

See [Optimization decision tree](@ref optimization-decision-tree).
"""
optimization_catalog() = copy(_OPTIMIZATION_CATALOG)

"""
    algorithm_impact(options)

Classify the selected run by the strongest scientific effect.  A physical
model change dominates a controlled numerical approximation, which dominates
roundoff-level computational rearrangements.

See [Evidence classes](@ref optimization-evidence-classes).
"""
function algorithm_impact(options::AlgorithmOptions)
    _check_algorithm_options(options)
    if options.retarded_real_part != :kramers_kronig ||
       options.self_energy_structure != :full ||
       options.transverse_momentum != :resolved ||
       options.embedding !== :full_cell_resolvent
        return :physical_model
    elseif options.localization in (:none, :real_space) ||
           options.hilbert === :product_integration ||
           options.energy_shift === :conservative_pair ||
           (
               options.contraction == :low_rank && (
                   options.low_rank_relative_tolerance > 0 ||
                   options.low_rank_maximum_rank > 0
               )
           ) ||
           options.kernel_build in (:tabulated, :adaptive_direct)
        return :controlled_numerical
    end
    return :physics_preserving
end

"""
    algorithm_manifest(options)

Return a flat dictionary suitable for YAML/HDF5 provenance and comparison
reports. Values use explicit strings instead of Julia implementation types.

See [YAML run configurations](https://github.com/Afonenko-QCL-NEGF/QCLNEGFRunner.jl/blob/main/docs/src/user/configuration.md) and
[Optimization decision tree](@ref optimization-decision-tree).
"""
function algorithm_manifest(options::AlgorithmOptions)
    _check_algorithm_options(options)
    return Dict{String,Any}(
        "solver_backend" => String(options.solver_backend),
        "energy_shift" => String(options.energy_shift),
        "hilbert" => String(options.hilbert),
        "contraction" => String(options.contraction),
        "kernel_build" => String(options.kernel_build),
        "mixing" => String(options.mixing),
        "localization" => String(options.localization),
        "embedding" => String(options.embedding),
        "embedding_periods" => options.embedding_periods,
        "occupation_normalization" => String(options.occupation_normalization),
        "retarded_real_part" => String(options.retarded_real_part),
        "self_energy_structure" => String(options.self_energy_structure),
        "transverse_momentum" => String(options.transverse_momentum),
        "low_rank_relative_tolerance" => options.low_rank_relative_tolerance,
        "low_rank_maximum_rank" => options.low_rank_maximum_rank,
        "anderson_history_depth" => options.anderson_history_depth,
        "anderson_damping" => options.anderson_damping,
        "anderson_regularization" => options.anderson_regularization,
        "impact" => String(algorithm_impact(options)),
    )
end
