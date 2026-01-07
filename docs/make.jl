using Documenter
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

makedocs(
    sitename = "BendersX.jl",
    format = Documenter.HTML(),
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
        "API" => "api.md",
    ],
    # sitelogo = "assets/logo.svg",
)