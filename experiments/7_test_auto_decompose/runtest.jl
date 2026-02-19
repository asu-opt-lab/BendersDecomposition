using BendersDecomposition
using Test
using JuMP


@testset "Auto Decompose Tests" begin
    @info "Running Auto Decompose Tests"
    # include("ufl_seq.jl")
    # include("ufl_in_out.jl")
    # include("ufl_callback.jl")
    # include("cfl_seq.jl")
    # include("cfl_in_out.jl")
    # include("cfl_callback.jl")
    include("MIPLIB.jl")
    @info "Auto Decompose Tests completed"
end
