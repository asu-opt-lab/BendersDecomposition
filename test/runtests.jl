using BendersDecomposition
using Test

@testset "BendersDecomposition.jl" begin
    # Write your tests here.
    data = UFLPData(
        5, 10, 
        rand(10) .* 50,     # Random demands
        rand(5) .* 100,     # Random fixed_costs
        rand(5, 10) .* 10,  # Random costs
    )
end
