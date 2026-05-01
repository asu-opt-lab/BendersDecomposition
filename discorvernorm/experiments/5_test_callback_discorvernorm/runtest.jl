using BendersX
using Test
using JuMP

isdefined(Main, :DiscorverNormDCGLP) || include(normpath(joinpath(@__DIR__, "..", "..", "src", "DiscorverNormDCGLP.jl")))
using .DiscorverNormDCGLP

@testset "Callback DiscorverNorm Tests" begin
    @info "Running DiscorverNormDCGLP callback tests"
    include("cfl.jl")
    @info "DiscorverNormDCGLP callback tests completed"
end
