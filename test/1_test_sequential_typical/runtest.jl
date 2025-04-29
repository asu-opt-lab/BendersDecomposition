using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition
import BendersDecomposition: generate_cuts


@testset "Sequential Typical Tests" begin
    @info "Running Sequential Typical Tests"
    include("ufl.jl")
    include("cfl.jl")
    include("scfl.jl")
    @info "Sequential Typical Tests completed"
end