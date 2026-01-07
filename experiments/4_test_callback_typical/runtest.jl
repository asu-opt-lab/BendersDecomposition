using BendersX
using Test
using JuMP

@testset "Callback Typical Tests" begin
    @info "Running Callback Typical Tests"
    include("ufl.jl")
    include("cfl.jl")
    include("scfl.jl")
    include("snip.jl")
    @info "Callback Typical Tests completed"
end