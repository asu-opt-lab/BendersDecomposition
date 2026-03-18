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
    format = Documenter.HTML(prettyurls=false, size_threshold_warn=150 * 2^10),
    checkdocs = :public,
    warnonly = [:cross_references],
    pages = [
        "Home" => "index.md",
        "Guides" => [
            "Getting Started" => "tutorials/getting_started.md",
            "Architecture" => "user_guide.md",
            "Modeling Guide" => "modeling_guide.md",
            "Problem Library" => "problem_library.md",
            "Oracle Guide" => "tutorials/oracles.md",
            "Environment Guide" => "tutorials/envs.md",
            "Extending BendersX" => "extending.md",
        ],
        "Reference" => [
            "CFLP Demo" => "tutorials/cflp_demo.md",
            "Experiments and Reproducibility" => "tutorials/examples.md",
            "API Reference" => "api.md",
        ],
    ],
)
