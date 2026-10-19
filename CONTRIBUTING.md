# Contributing

For documentation-only text or link corrections, review the diff and check the
changed links, anchors and formatting. This does not require Julia, package
installation, a full test suite or documentation generation. Record numerical and
rendering checks that were not run. Executable examples, documentation build
recipes and public API changes need the relevant checks for their behavior.

For code changes, run relevant targeted checks in the prepared environment.
`deno task check` checks the TypeScript orchestration; it does not test numerical
behavior. `deno task test` runs Julia tests, and `deno task docs` runs the Julia
documentation generator; both may instantiate dependencies when no prepared
project is selected. Environment preparation, the full suite and documentation
generation belong to their assigned stages with a finite shared budget. Missing
runtime is a reported check limitation, not permission to mark a check as passed.

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

Integration CI is defined in the [qcl-negf superproject](https://github.com/Afonenko-QCL-NEGF/qcl-negf) and uses its local runner. Update the component gitlink there to check a change with the complete selected source graph.
