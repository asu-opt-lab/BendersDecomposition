using Test
using BendersX

@testset "BendersX.jl" begin
    include("test_scalarize_constraints.jl")
    include("test_unified_oracle.jl")
    include("test_pareto_oracle.jl")
end