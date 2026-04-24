using BendersX
using Test
using JuMP

isdefined(Main, :DirectionalPolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "DirectionalPolarDCGLP.jl")))
using .DirectionalPolarDCGLP

@testset "Callback Directional Disjunctive Tests" begin
    @info "Running DirectionalPolarDCGLP callback disjunctive tests"
    include("cfl.jl")
    @info "DirectionalPolarDCGLP callback disjunctive tests completed"
end
