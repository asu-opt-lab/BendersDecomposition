using Documenter
using BendersDecomposition
using Literate

const _TUTORIAL_DIR = joinpath(@__DIR__, "src", "tutorial")
const _EXAMPLES_DIR = joinpath(@__DIR__, "src", "examples")

function literate_directory(dir)
    for filename in filter(f -> endswith(f, ".jl"), readdir(dir))
        Literate.markdown(joinpath(dir, filename), dir; documenter = true)
    end
end

literate_directory(_TUTORIAL_DIR)
literate_directory(_EXAMPLES_DIR)

makedocs(
    sitename = "BendersDecomposition.jl",
    format = Documenter.HTML(),
    modules = [BendersDecomposition],
    authors = "Kaiwen Fang",
    pages = [
        "Home" => "index.md",
        "Tutorials" => [
            "tutorial/introduction.md",
            "tutorial/beginner.md",
            "tutorial/advanced.md",
        ],
        "Examples" => [
            "examples/intro.md",
            "examples/CFLP.md",
            "examples/UFLP.md",
            "examples/MCNDP.md",
            "examples/SCFLP.md",
            "examples/SNIP.md",
        ],
        "API Reference" => "api.md"
    ]
)



