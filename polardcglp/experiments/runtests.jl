using BendersX
using Test

isdefined(Main, :SimplexNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "src", "SimplexNormDCGLP.jl")))
using .SimplexNormDCGLP

@testset "SimplexNormDCGLP Experiments" begin
    include("3_test_sequential_simplex_norm_dcglp/runtest.jl")
    include("3_test_sequential_directional_polar/runtest.jl")
    include("3_test_sequential_vertical_reverse_polar/runtest.jl")
    include("3_test_sequential_simplex_norm_test/runtest.jl")
end
