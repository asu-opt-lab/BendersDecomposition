using BendersX
using Test
using JuMP

isdefined(Main, :PolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "PolarDCGLP.jl")))
using .PolarDCGLP

@testset "Callback Directional Disjunctive Tests" begin
    @info "Running DirectionalPolarDCGLP callback disjunctive tests"
    include("cfl.jl")
    @info "DirectionalPolarDCGLP callback disjunctive tests completed"
end
