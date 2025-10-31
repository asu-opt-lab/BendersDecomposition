using BendersDecomposition
using Test
using JuMP

@testset "Callback Disjunctive Tests" begin
    @info "Running Callback Disjunctive Tests"
    include("ufl.jl")
    include("cfl.jl")
    include("scfl.jl")
    @info "Callback Disjunctive Tests completed"
end