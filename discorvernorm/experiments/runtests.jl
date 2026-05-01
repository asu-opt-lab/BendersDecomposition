using BendersX
using Test

isdefined(Main, :DiscorverNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "src", "DiscorverNormDCGLP.jl")))
using .DiscorverNormDCGLP

@testset "DiscorverNorm Experiments" begin
    include("1_test_sequential_discorvernorm/runtest.jl")
    include("2_test_sequential_in_out_discorvernorm/runtest.jl")
    include("3_test_sequential_discorvernorm/runtest.jl")
    include("5_test_callback_discorvernorm/runtest.jl")
end
