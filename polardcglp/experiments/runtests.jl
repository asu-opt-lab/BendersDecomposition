using BendersX
using Test

isdefined(Main, :PolarDCGLP) || include(normpath(joinpath(@__DIR__, "..", "src", "PolarDCGLP.jl")))
using .PolarDCGLP

@testset "PolarDCGLP Experiments" begin
    # include("3_test_sequential_disjunctive/runtest.jl")
    # include("5_test_callback_disjunctive/runtest.jl")
    # include("3_test_sequential_directional_disjunctive/runtest.jl")
    include("5_test_callback_directional_disjunctive/runtest.jl")
end
