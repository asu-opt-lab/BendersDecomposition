using BendersX
using Test
using JuMP

isdefined(Main, :DiscorverNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "DiscorverNormDCGLP.jl")))
using .DiscorverNormDCGLP

@testset "Sequential DiscorverNorm Tests" begin
    @info "Running DiscorverNormDCGLP sequential tests"
    include("cfl.jl")
    @info "DiscorverNormDCGLP sequential tests completed"
end
