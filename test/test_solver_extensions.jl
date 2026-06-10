using BendersX
using Pkg
using Test

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const PROJECT_TOML = Pkg.TOML.parsefile(joinpath(REPO_ROOT, "Project.toml"))

@testset "Optional Solver Extensions" begin
    @test haskey(PROJECT_TOML, "weakdeps")
    @test haskey(PROJECT_TOML, "extensions")
    @test get(PROJECT_TOML["extensions"], "BendersXCPLEXExt", nothing) == "CPLEX"
    @test haskey(PROJECT_TOML["weakdeps"], "CPLEX")

    @test isfile(joinpath(REPO_ROOT, "ext", "BendersXCPLEXExt.jl"))
    @test Base.get_extension(BendersX, :BendersXCPLEXExt) === nothing
end
