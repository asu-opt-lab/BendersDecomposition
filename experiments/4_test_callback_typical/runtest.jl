using BendersDecomposition
using Test
using JuMP

@testset "Callback Typical Tests" begin
    @info "Running Callback Typical Tests"
    include("ufl.jl")
    include("cfl.jl")
    include("scfl.jl")
    @info "Callback Typical Tests completed"
end