using BendersDecomposition
using Test
using JuMP

@testset "Sequential Typical Tests" begin
    @info "Running Sequential InOut Typical Tests"
    include("ufl.jl")
    include("cfl.jl")
    include("scfl.jl")
    # include("snip.jl")
    @info "Sequential InOut Typical Tests completed"
end