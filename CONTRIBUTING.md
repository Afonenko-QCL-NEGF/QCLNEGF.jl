# Contributing

Run `deno task check`, `deno task test` and `deno task docs` before proposing a change.
Keep the public model dimensional and document normalization, index order and
approximation assumptions. A numerical change needs an independent analytic or
physical oracle; agreement between two implementations is supplementary evidence.
Do not loosen acceptance thresholds to make a failed calculation pass.

Keep reference operators independent from optimized operators. Preserve the
separation between convergence, discretization evidence and experimental validation.
Changes to public types or solver callbacks must be tested with QCLNEGFRunner.jl.
Filesystem access, YAML, persistence, resource discovery and presentation belong to
that downstream package. Core numerical routines accept explicit inputs and return
in-memory results.

Integration CI is defined in the [qcl-negf superproject](https://github.com/AfonenkoA/qcl-negf) and uses its local runner. Update the component gitlink there to check a change with the complete selected source graph.
