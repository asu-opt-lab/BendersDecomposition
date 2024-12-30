using Test
using JuMP
# using Gurobi

# Import the necessary modules and functions
using BendersDecomposition

@testset "SCFLP Master Problem" begin
    # Create a mock data structure for testing
    # Create mock data
    n_facilities = 3
    n_customers = 4
    n_scenarios = 2
    ratio = 2

    # Generate test data
    data = generate_stochastic_capacited_facility_location(
        n_facilities,
        n_customers,
        n_scenarios,
        ratio
    )

    # Test StandardSCFLPMasterProblem
    @testset "StandardSCFLPMasterProblem" begin
        cut_strategy = ClassicalCut()
        mp = create_master_problem(data, cut_strategy)

        @test mp isa SCFLPMasterProblem
        @test mp.model isa Model
        @test haskey(mp.var, :x)
        @test haskey(mp.var, :t)
        @test length(mp.var[:x]) == data.n_facilities
        @test length(mp.var[:t]) == data.n_scenarios
        @test mp.obj_value == 0.0

        # Check if the objective function is set correctly
        obj = objective_function(mp.model)
        @test obj isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        @test length(obj.terms) == data.n_facilities + data.n_scenarios  # x variables + t variable
    end

end
