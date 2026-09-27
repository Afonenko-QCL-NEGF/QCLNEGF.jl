# [Architecture](@id developer-architecture)

QCLNEGF owns the in-memory scientific model. Dependencies point from
[QCLNEGFRunner.jl](https://github.com/Afonenko-QCL-NEGF/QCLNEGFRunner.jl) into QCLNEGF;
this package does not import its runner.

## [Source ownership](@id developer-source-map)

| Source | Responsibility |
|---|---|
| `src/domain` | Dimensional parameters, physical models and convergence contracts |
| `src/numerics/reference` | Independent literal operators and their algebraic definitions |
| `src/numerics/optimized` | Optimized kernels, bounded workspaces and nonlinear solve |
| `src/physics` | Physical observables and optical-response operators |
| `src/application` | Numerical algorithm descriptors and model-capability metadata |
| `src/composition` | Operating-point reuse and numerical precompilation workload |

A caller constructs `PhysicalParameters`, `NumericalParameters`, scattering and
solver options, builds a problem, and receives an in-memory solution. Explicit
callbacks can observe or interrupt the solve without choosing a storage format.
Thread width and algorithm options are explicit numerical inputs. Resource discovery
and process admission are Runner responsibilities.

YAML, HDF5, checkpoints, result directories, reports and CairoMakie integration are
owned by Runner. AiiDA owns durable distributed workflow state; Slurm owns CPU and
memory allocation. The mathematical acceptance criteria are independent of those
execution choices.

See [operator contracts](@ref physics-operators), [traceability](@ref traceability),
and [convergence policy](@ref convergence-policy).
