using BendersX
using Test
using JuMP

@testset "Sequential Typical FLCAP Tests" begin
    @info "Running sequential typical FLCAP tests"
    include("cfl.jl")
    @info "Sequential typical FLCAP tests completed"
end
