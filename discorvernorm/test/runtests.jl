using Test
using BendersX

isdefined(Main, :DiscorverNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "src", "DiscorverNormDCGLP.jl")))
using .DiscorverNormDCGLP

include("test_discorvernorm_dcglp.jl")
