using BendersDecomposition
using Test
using JuMP



@testset "Sequential Typical Tests" begin
    @info "Running Sequential Typical Tests"
    include("ufl.jl")
    include("cfl.jl")
    include("scfl.jl")
    @info "Sequential Typical Tests completed"
end