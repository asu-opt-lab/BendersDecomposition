using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition
import BendersDecomposition: generate_cuts

@testset "Callback Typical Tests" begin
    @info "Running Callback Typical Tests"
    include("ufl.jl")
    include("cfl.jl")
    include("scfl.jl")
    @info "Callback Typical Tests completed"
end