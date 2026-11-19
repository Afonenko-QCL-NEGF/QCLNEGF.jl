using QCLNEGF, Documenter
DocMeta.setdocmeta!(QCLNEGF, :DocTestSetup, :(using QCLNEGF, Unitful, LinearAlgebra); recursive = true)
makedocs(
    modules = [QCLNEGF],
    remotes = nothing,
    sitename = "QCLNEGF",
    format = Documenter.HTML(prettyurls = false, edit_link = nothing, lang = "ru",
        footer = "Собрано с помощью [Documenter.jl](https://documenter.juliadocs.org/) и [Julia](https://julialang.org/).",
        repolink = "https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl"),
    pages = [
        "Главная" => "index.md",
        "Использование" => ["user/convergence.md"],
        "API" => ["api/operators.md", "api/public.md", "api/extensions.md"],
        "Теория" => ["theory/01_model.md", "theory/02_reference2019.md", "theory/03_arrays.md", "theory/04_scaling.md", "theory/05_basis.md", "theory/06_periodicity.md", "theory/07_greens.md", "theory/08_kernels.md", "theory/09_scba.md", "theory/10_poisson.md", "theory/11_observables.md", "theory/12_validation.md", "theory/13_visualization.md", "theory/14_cost.md", "theory/15_limitations.md", "theory/16_input_contract.md", "theory/17_implementation_plan.md", "theory/18_optical_response.md", "theory/19_production.md", "theory/20_optimization_decision_tree.md"],
        "Модели и численные критерии" => ["physics/models.md", "physics/numerical-policy.md"],
        "Физические иллюстрации" => ["tutorials/localization.md", "tutorials/physical-atlas.md"],
        "Разработка" => ["developer/architecture.md", "developer/physics_operators.md", "developer/traceability.md"],
        "Справочные материалы" => ["glossary.md", "references.md", "references_optimization.md"],
    ],
    doctest = true,
    checkdocs = :exports,
    warnonly = false,
)

include("localize_html.jl")
localize_documenter!(joinpath(@__DIR__, "build"))
