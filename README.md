# QCLNEGF.jl

Julia **1.13.0** library for nonequilibrium Green-function transport in layered
quantum-cascade heterostructures. Version **0.2.0** contains dimensional model
types, independent reference and optimized numerical operators, SCBA and Poisson
solvers, physical observables, and optical-response calculations.

The library solves an explicitly constructed problem in memory. YAML studies,
files, checkpoints, reports, plots and command-line execution belong to
[QCLNEGFRunner.jl](https://github.com/AfonenkoA/QCLNEGFRunner.jl).

## Use

From a checkout, instantiate the Julia project:

```console
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

A small direct calculation uses a declared educational preset:

```julia
using QCLNEGF

problem = build_problem(
    physical = reference_parameters(),
    numerical = tutorial_numerics(),
    scattering = ScatteringOptions(LO=false, acoustic=true,
        impurity=false, IFR=false, alloy=false),
)
solution = solve(problem; options=tutorial_options())
solution.status
```

Inspect solution status, residuals and physical gates before interpreting an
observable. A completed iteration sequence is not evidence of grid convergence or
experimental agreement. The small tutorial grid demonstrates the equations; it is
not a quantitative device prediction.

The reference geometry comes from “Thermoelectrically cooled THz quantum cascade
laser operating up to 210 K”, [DOI 10.1063/1.5110305](https://doi.org/10.1063/1.5110305).
General layer sequences are supported through `PhysicalParameters`.

## Development

```console
deno task check
deno task test
deno task docs
```

The checked dependency manifest and exact Julia version define the development
environment. Numerical tests use bounds checking, independent analytic cases,
conservation identities and comparisons between literal and optimized operators.
CI runs on trusted local GitHub Actions runners managed by
[qcl-negf-platform](https://github.com/AfonenkoA/qcl-negf-platform).

Start with the [model](docs/src/theory/01_model.md),
[in-memory workflow](docs/src/theory/17_implementation_plan.md),
[operator reference](docs/src/api/operators.md), and
[architecture](docs/src/developer/architecture.md). Mathematical reference chapters
are in Russian. Build the browsable documentation with `deno task docs`.

[QCLNEGFRunner.jl](https://github.com/AfonenkoA/QCLNEGFRunner.jl) consumes this API.
[qcl-negf-research](https://github.com/AfonenkoA/qcl-negf-research) owns study definitions;
[qcl-negf-aiida](https://github.com/AfonenkoA/qcl-negf-aiida) and Slurm schedule them.

License: MIT.
