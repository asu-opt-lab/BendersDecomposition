using Documenter
using BendersX

makedocs(
    modules = [BendersX],
    sitename = "BendersX.jl",
    format = Documenter.HTML(prettyurls=false),
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
