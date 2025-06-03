using Test
using JuMP
using CPLEX
using Printf
using DataFrames
using Logging
using BendersDecomposition
import BendersDecomposition: generate_cuts


@testset "Callback Disjunctive Tests" begin
    @info "Running Callback Disjunctive Tests"
    include("ufl.jl")
    include("cfl.jl")
    include("scfl.jl")
    include("snip.jl")
    @info "Callback Disjunctive Tests completed"
end