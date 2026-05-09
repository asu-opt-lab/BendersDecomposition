using BendersX
using Test
using JuMP

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

@testset "Sequential Vertical Reverse Polar Tests" begin
    @info "Running VerticalReversePolarDCGLP sequential disjunctive tests"
    include("ufl.jl")
    include("cfl.jl")
    @info "VerticalReversePolarDCGLP sequential disjunctive tests completed"
end
