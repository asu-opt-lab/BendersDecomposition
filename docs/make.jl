using Documenter
using BendersX

makedocs(
    sitename = "BendersX.jl",
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        # "Tutorial" => [
        #     "Data" => "Tutorial/data.md",
        #     "Modeling Interface" => "Tutorial/modeling.md",
        #     "Running Benders" => "Tutorial/envs.md",
        #     "Example - Capacitated Facility Location Problem (CFLP)" => "Tutorial/example.md",
        # ],
        # "API" => "api.md",
    ]
)