using BendersX
using Test

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

@testset "SimplexNormDCGLP Experiments" begin
    # include("3_test_sequential_disjunctive/runtest.jl")
    # include("5_test_callback_disjunctive/runtest.jl")
    # include("3_test_sequential_directional_disjunctive/runtest.jl")
    include("5_test_callback_directional_disjunctive/runtest.jl")
end
