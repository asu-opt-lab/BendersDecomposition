using Test
using JuMP
# using Gurobi

# Import the necessary modules and functions
using BendersDecomposition

@testset "CFLP Master Problem" begin
    # Create a mock data structure for testing
    # Create mock data
    data = CFLPData(
        5, 10, 
        rand(5) .* 200,     # Random capacities
        rand(10) .* 50,     # Random demands
        rand(5) .* 100,     # Random fixed_costs
        rand(5, 10) .* 10,  # Random costs
    )
    # Test StandardCFLPMasterProblem
    @testset "StandardCFLPMasterProblem" begin
        cut_strategy = ClassicalCut()
        mp = create_master_problem(data, cut_strategy)

        @test mp isa CFLPMasterProblem
        @test mp.model isa Model
        @test haskey(mp.var, :x)
        @test haskey(mp.var, :t)
        @test length(mp.var[:x]) == data.n_facilities
        @test length(mp.var[:t]) == 1
        @test mp.obj_value == 0.0
        @test mp.t_value == 0.0
        @test mp.x_value == zeros(data.n_facilities)
        
        # Check if the objective function is set correctly
        obj = objective_function(mp.model)
        @test obj isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        @test length(obj.terms) == data.n_facilities + 1  # x variables + t variable
    end

end
