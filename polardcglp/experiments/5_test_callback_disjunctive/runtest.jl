using BendersX
using Test
using JuMP

isdefined(Main, :PolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "PolarDCGLP.jl")))
using .PolarDCGLP

@testset "Callback Disjunctive Tests" begin
    @info "Running PolarDCGLP callback disjunctive tests"
    include("ufl.jl")
    include("cfl.jl")
    @info "PolarDCGLP callback disjunctive tests completed"
end
