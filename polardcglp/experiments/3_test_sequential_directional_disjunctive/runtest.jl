using BendersX
using Test
using JuMP

isdefined(Main, :PolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "PolarDCGLP.jl")))
using .PolarDCGLP

@testset "Sequential Directional Disjunctive Tests" begin
    @info "Running DirectionalPolarDCGLP sequential disjunctive tests"
    include("cfl.jl")
    @info "DirectionalPolarDCGLP sequential disjunctive tests completed"
end
