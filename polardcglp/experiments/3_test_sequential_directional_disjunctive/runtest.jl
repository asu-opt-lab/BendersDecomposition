using BendersX
using Test
using JuMP

isdefined(Main, :DirectionalPolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "DirectionalPolarDCGLP.jl")))
using .DirectionalPolarDCGLP

@testset "Sequential Directional Disjunctive Tests" begin
    @info "Running DirectionalPolarDCGLP sequential disjunctive tests"
    include("cfl.jl")
    @info "DirectionalPolarDCGLP sequential disjunctive tests completed"
end
