using BendersX
using Test
using JuMP

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

@testset "Sequential Directional Disjunctive Tests" begin
    @info "Running DirectionalPolarDCGLP sequential disjunctive tests"
    include("cfl.jl")
    @info "DirectionalPolarDCGLP sequential disjunctive tests completed"
end
