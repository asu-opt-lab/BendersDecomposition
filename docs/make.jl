using Pkg

Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()

using BendersX
using Documenter
using Literate

const TUTORIALS_SRC = joinpath(@__DIR__, "src", "tutorials")
const TUTORIALS_OUT = joinpath(@__DIR__, "src", "tutorials")

const LITERATE_FILES = [
    "cflp_demo.jl",
]

for file in LITERATE_FILES
    input_file = joinpath(TUTORIALS_SRC, file)
    if isfile(input_file)
        println("Processing Literate file: ", file)
        Literate.markdown(
            input_file,
            TUTORIALS_OUT;
            documenter = true,
            execute = false,
        )
    else
        @warn "Literate source file not found: $input_file"
    end
end

makedocs(
    modules = [BendersX],
    sitename = "BendersX.jl",
    format = Documenter.HTML(prettyurls=false),
    checkdocs = :exports,
    warnonly = [:cross_references],
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "Getting Started" => "tutorials/getting_started.md",
            "CFLP Demo" => "tutorials/cflp_demo.md",
            "Swapping Oracles and Adjusting Their Behaviors" => "tutorials/oracles.md",
            "Swapping Environments and Adjusting Their Behaviors" => "tutorials/envs.md",
            "Examples" => "tutorials/examples.md"
        ],
        "User Guide" => "user_guide.md",
        "API" => "api.md"
    ],
)
