# Ordinary offline documentation uses a fresh Root-owned attempt namespace.
# These local guards validate layout and inputs before imports; Root owns staged
# payload/assets/cache closure. pilot-v1 remains a compatibility profile.
_qcl_docs_profile = get(ENV, "QCL_DOCS_PROFILE", "")
if !isempty(_qcl_docs_profile) && !(_qcl_docs_profile in ("offline-v1", "pilot-v1"))
    error("Unsupported QCL_DOCS_PROFILE")
end

if _qcl_docs_profile in ("offline-v1", "pilot-v1")
    # Base only before all guards complete: no imports, setup or output mutation.
    function _qcl_docs_path(name)
        value = get(ENV, name, "")
        isempty(value) && error("Missing required documentation path: " * name)
        startswith(value, "/") && isabspath(value) ||
            error("Documentation path must be Linux absolute: " * name)
        (endswith(value, "/") || occursin("//", value) ||
         any(part -> part == "." || part == "..", split(value, '/'))) &&
            error("Documentation path has noncanonical components: " * name)
        value == normpath(value) || error("Documentation path must be normalized: " * name)
        cursor = value
        while true
            islink(cursor) && error("Documentation path has a symlink ancestor: " * name)
            parent = dirname(cursor)
            parent == cursor && break
            cursor = parent
        end
        return value
    end

    function _qcl_docs_directory(path)
        isdir(path) || error("Required documentation directory is missing")
        realpath(path) == path || error("Documentation directory is not canonical")
        return path
    end

    function _qcl_docs_contains(parent, child)
        parent_parts = splitpath(parent)
        child_parts = splitpath(child)
        return length(child_parts) >= length(parent_parts) &&
            child_parts[1:length(parent_parts)] == parent_parts
    end

    _qcl_docs_root = _qcl_docs_path("QCL_DOCS_ROOT")
    _qcl_docs_attempt = _qcl_docs_path("QCL_DOCS_ATTEMPT")
    _qcl_docs_source = _qcl_docs_path("QCL_DOCS_SOURCE")
    _qcl_docs_build = _qcl_docs_path("QCL_DOCS_BUILD")
    _qcl_docs_directory(_qcl_docs_root)
    _qcl_docs_directory(_qcl_docs_attempt)
    _qcl_docs_directory(_qcl_docs_source)
    _qcl_docs_directory(dirname(_qcl_docs_build))
    (_qcl_docs_contains(_qcl_docs_root, _qcl_docs_attempt) ||
     _qcl_docs_contains(_qcl_docs_attempt, _qcl_docs_root)) &&
        error("Documentation ROOT and ATTEMPT must not overlap")
    (_qcl_docs_contains(_qcl_docs_source, _qcl_docs_build) ||
     _qcl_docs_contains(_qcl_docs_build, _qcl_docs_source)) &&
        error("Documentation SOURCE and BUILD must not overlap")
    (ispath(_qcl_docs_build) || islink(_qcl_docs_build)) &&
        error("Documentation BUILD must not already exist; use a fresh attempt")
    VERSION == v"1.13.0" || error("Documentation requires Julia 1.13.0")
    _qcl_docs_project = joinpath(_qcl_docs_root, "julia", "Project.toml")
    _qcl_docs_manifest = joinpath(_qcl_docs_root, "julia", "Manifest.toml")
    # Check all project ancestors and leaf files as explicit paths, too.
    for (name, expected) in (("Project.toml", _qcl_docs_project),
                             ("Manifest.toml", _qcl_docs_manifest))
        cursor = expected
        while true
            islink(cursor) && error("Documentation project path has a symlink ancestor: " * name)
            parent = dirname(cursor)
            parent == cursor && break
            cursor = parent
        end
        isfile(expected) || error("Missing regular documentation project input: " * name)
        realpath(expected) == expected || error("Documentation project input is not canonical")
    end
    Base.active_project() == _qcl_docs_project ||
        error("Documentation requires ROOT/julia/Project.toml as active project")
    _qcl_docs_source == joinpath(_qcl_docs_attempt, "stage", "core", "src") ||
        error("Documentation SOURCE does not match the fixed core layout")
    _qcl_docs_build == joinpath(_qcl_docs_attempt, "raw", "core") ||
        error("Documentation BUILD does not match the fixed core layout")
    joinpath(realpath(dirname(_qcl_docs_build)), "core") == _qcl_docs_build ||
        error("Documentation BUILD candidate is not canonical")
    # Root separately rejects symlink/special-file payload entries and binds all
    # staged bytes, effect manifests and assets/licenses for the ordinary build.

    using QCLNEGF, Documenter
    makedocs(
        root = _qcl_docs_root,
        source = _qcl_docs_source,
        build = _qcl_docs_build,
        modules = [QCLNEGF],
        remotes = nothing,
        sitename = "QCLNEGF",
        format = Documenter.HTML(prettyurls = false, edit_link = nothing,
            disable_git = true, analytics = "", warn_outdated = false,
            # Совместимо с Documenter 1.17; Root заменяет raw JS и рендерит формулы локально.
            prerender = false, mathengine = Documenter.KaTeX(),
            assets = ["assets/qcl-fixed-light.css"],
            repolink = "https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl"),
        pages = [
            "Home" => "index.md",
            "Using the package" => ["user/convergence.md"],
            "API" => ["api/operators.md", "api/public.md"],
            "Theory" => ["theory/01_model.md", "theory/02_reference2019.md", "theory/03_arrays.md", "theory/04_scaling.md", "theory/05_basis.md", "theory/06_periodicity.md", "theory/07_greens.md", "theory/08_kernels.md", "theory/09_scba.md", "theory/10_poisson.md", "theory/11_observables.md", "theory/12_validation.md", "theory/13_visualization.md", "theory/14_cost.md", "theory/15_limitations.md", "theory/16_input_contract.md", "theory/17_implementation_plan.md", "theory/18_optical_response.md", "theory/19_production.md", "theory/20_optimization_decision_tree.md"],
            "Models and policy" => ["physics/models.md", "physics/numerical-policy.md"],
            "Physical illustrations" => ["tutorials/localization.md", "tutorials/physical-atlas.md"],
            "Development" => ["developer/architecture.md", "developer/physics_operators.md", "developer/traceability.md"],
            "Reference material" => ["glossary.md", "references.md", "references_optimization.md"],
        ],
        pagesonly = false,
        doctest = false,
        linkcheck = false,
        checkdocs = :exports,
        warnonly = false,
        draft = false,
    )
else
using QCLNEGF, Documenter
DocMeta.setdocmeta!(QCLNEGF, :DocTestSetup, :(using QCLNEGF, Unitful, LinearAlgebra); recursive = true)
makedocs(
    modules = [QCLNEGF],
    remotes = nothing,
    sitename = "QCLNEGF",
    format = Documenter.HTML(prettyurls = false, edit_link = nothing,
        repolink = "https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl"),
    pages = [
        "Home" => "index.md",
        "Using the package" => ["user/convergence.md"],
        "API" => ["api/operators.md", "api/public.md"],
        "Theory" => ["theory/01_model.md", "theory/02_reference2019.md", "theory/03_arrays.md", "theory/04_scaling.md", "theory/05_basis.md", "theory/06_periodicity.md", "theory/07_greens.md", "theory/08_kernels.md", "theory/09_scba.md", "theory/10_poisson.md", "theory/11_observables.md", "theory/12_validation.md", "theory/13_visualization.md", "theory/14_cost.md", "theory/15_limitations.md", "theory/16_input_contract.md", "theory/17_implementation_plan.md", "theory/18_optical_response.md", "theory/19_production.md", "theory/20_optimization_decision_tree.md"],
        "Models and policy" => ["physics/models.md", "physics/numerical-policy.md"],
        "Physical illustrations" => ["tutorials/localization.md", "tutorials/physical-atlas.md"],
        "Development" => ["developer/architecture.md", "developer/physics_operators.md", "developer/traceability.md"],
        "Reference material" => ["glossary.md", "references.md", "references_optimization.md"],
    ],
    doctest = true,
    checkdocs = :exports,
    warnonly = false,
)
end
