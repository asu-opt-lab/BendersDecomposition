using Documenter
using Literate
using BendersX

println("Module loaded. Exports:", names(BendersX, all=false))
for n in names(BendersX, all=true)
    if isdefined(BendersX, n)
        obj = getfield(BendersX, n)
        has = Base.Docs.doc(obj) !== nothing
        if has
            println("DOC: ", n)
        end
    end
end

# Process Literate.jl source files
const TUTORIALS_SRC = joinpath(@__DIR__, "src", "tutorials")
const TUTORIALS_OUT = joinpath(@__DIR__, "src", "tutorials")

# List of Literate.jl source files to process
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
            execute = false,  # Set to true if you want code execution
        )
    else
        @warn "Literate source file not found: $input_file"
    end
end

makedocs(
    sitename = "BendersX.jl",
    format = Documenter.HTML(prettyurls=false),
    pages = [
        "Home" => "index.md",
        # "Tutorials" => "tutorial.md",
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
    # sitelogo = "assets/logo.svg",
)