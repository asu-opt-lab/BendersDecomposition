using Documenter
using BendersX

makedocs(
    modules = [BendersX],
    sitename = "BendersX.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://asu-opt-lab.github.io/BendersDecomposition/",
    ),
    checkdocs = :public,
    pages = [
        "Home" => "index.md",
        # "Tutorials" => "tutorial.md",
        "Tutorials" => [
            "Getting Started" => "tutorials/getting_started.md",
            "Swapping Oracles and Adjusting Their Behaviors" => "tutorials/oracles.md",
            "Swapping Environments and Adjusting Their Behaviors" => "tutorials/envs.md",
            "Examples" => "tutorials/examples.md"
        ],
        "User Guide" => "user_guide.md",
        "API" => "api.md"
    ],
    # sitelogo = "assets/logo.svg",
)

deploydocs(
    repo = "github.com/asu-opt-lab/BendersDecomposition.git",
    devbranch = "main",
    versions = nothing,
)
