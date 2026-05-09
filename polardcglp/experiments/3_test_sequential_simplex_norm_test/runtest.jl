using BendersX
using Test
using JuMP

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

@testset "Sequential Simplex Norm Test Tests" begin
    @info "Running SimplexNormTestDCGLP sequential disjunctive tests"
    include("ufl.jl")
    include("cfl.jl")
    @info "SimplexNormTestDCGLP sequential disjunctive tests completed"
end
