using BendersX
using Test
using JuMP

@testset "Sequential InOut Typical FLCAP Tests" begin
    @info "Running sequential InOut typical FLCAP tests"
    include("cfl.jl")
    @info "Sequential InOut typical FLCAP tests completed"
end
