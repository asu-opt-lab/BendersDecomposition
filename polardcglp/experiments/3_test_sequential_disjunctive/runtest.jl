using BendersX
using Test
using JuMP

isdefined(Main, :PolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "PolarDCGLP.jl")))
using .PolarDCGLP

@testset "Sequential Disjunctive Tests" begin
    @info "Running PolarDCGLP sequential disjunctive tests"
    include("ufl.jl")
    include("cfl.jl")
    @info "PolarDCGLP sequential disjunctive tests completed"
end
