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
        
        # Test initial values of x_value and t_value
        @test length(mp.x_value) == data.n_facilities
        @test all(mp.x_value .== 0.0)
        @test length(mp.t_value) == data.n_scenarios
        @test all(mp.t_value .== 0.0)

        # Check if the objective function is set correctly
        obj = objective_function(mp.model)
        @test obj isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        @test length(obj.terms) == data.n_facilities + data.n_scenarios

        # Test t variable bounds
        for t_var in mp.var[:t]
            @test has_lower_bound(t_var)
            @test lower_bound(t_var) == -1e6
        end

        # Test capacity constraint
        max_demand = maximum(sum(demands) for demands in data.demands)
        cons = all_constraints(mp.model, AffExpr, MOI.GreaterThan{Float64})
        @test length(cons) == 1  # Should have one capacity constraint
        @test normalized_coefficient(cons[1], mp.var[:x][1]) == data.capacities[1]
    end

end
