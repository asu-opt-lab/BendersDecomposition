using BendersX
using Test
using JuMP

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

@testset "Callback Disjunctive Tests" begin
    @info "Running SimplexNormDCGLP callback disjunctive tests"
    include("ufl.jl")
    include("cfl.jl")
    @info "SimplexNormDCGLP callback disjunctive tests completed"
end
