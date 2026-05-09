using BendersX
using Test
using JuMP

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

@testset "Callback Directional Disjunctive Tests" begin
    @info "Running DirectionalPolarDCGLP callback disjunctive tests"
    include("cfl.jl")
    @info "DirectionalPolarDCGLP callback disjunctive tests completed"
end
