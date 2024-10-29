using Test
using JuMP
using Gurobi

# Import the necessary modules and functions
using BendersDecomposition

@testset "UFLP Master Problem" begin
    # Create a mock data structure for testing
    # Create mock data
    data = UFLPData(
        5, 10, 
        rand(10) .* 50,     # Random demands
        rand(5) .* 100,     # Random fixed_costs
        rand(5, 10) .* 10,  # Random costs
    )
    # Test StandardUFLPMasterProblem
    @testset "StandardUFLPMasterProblem" begin
        cut_strategy = StandardCut()
        mp = create_master_problem( data, cut_strategy)

        @test mp isa UFLPMasterProblem
        @test mp.model isa Model
        @test haskey(mp.var, :x)
        @test haskey(mp.var, :t)
        @test length(mp.var[:x]) == data.n_facilities
        @test mp.obj_value == 0.0

        # Check if the objective function is set correctly
        obj = objective_function(mp.model)
        @test obj isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        @test length(obj.terms) == data.n_facilities + 1  # x variables + t variable
    end

    # Test KnapsackUFLPMasterProblem
    @testset "FatKnapsackUFLPMasterProblem" begin
        cut_strategy = FatKnapsackCut()
        mp = create_master_problem(data, cut_strategy)

        @test mp isa UFLPMasterProblem
        @test mp.model isa Model
        @test haskey(mp.var, :x)
        @test haskey(mp.var, :t)
        @test length(mp.var[:x]) == data.n_facilities
        @test length(mp.var[:t]) == data.n_customers
        @test mp.obj_value == 0.0

        # Check if the objective function is set correctly
        obj = objective_function(mp.model)
        @test obj isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        @test length(obj.terms) == data.n_facilities + data.n_customers
    end

    @testset "SlimKnapsackUFLPMasterProblem" begin
        cut_strategy = SlimKnapsackCut()
        mp = create_master_problem(data, cut_strategy)

        @test mp isa UFLPMasterProblem
        @test mp.model isa Model
        @test haskey(mp.var, :x)
        @test haskey(mp.var, :t)
        @test length(mp.var[:x]) == data.n_facilities
        @test length(mp.var[:t]) == 1
        @test mp.obj_value == 0.0
    end
end
