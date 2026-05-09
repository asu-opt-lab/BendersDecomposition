using BendersX
using Test
using JuMP

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

@testset "Sequential Disjunctive Tests" begin
    @info "Running SimplexNormDCGLP sequential disjunctive tests"
    include("ufl.jl")
    include("cfl.jl")
    @info "SimplexNormDCGLP sequential disjunctive tests completed"
end
