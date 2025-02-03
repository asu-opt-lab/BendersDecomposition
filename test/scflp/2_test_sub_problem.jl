using Test
using JuMP
using Gurobi
using BendersDecomposition

@testset "SCFLP Sub Problem" begin
    data = read_stochastic_capacited_facility_location_problem("f10-c10-s5-r3-1")

    @testset "create_scenario_data" begin
        scenario_data = BendersDecomposition.create_scenario_data(data, 1)
        @test scenario_data isa CFLPData
        @test scenario_data.n_facilities == data.n_facilities
        @test scenario_data.n_customers == data.n_customers
        @test scenario_data.capacities == data.capacities
        @test scenario_data.demands == data.demands[1]
        @test scenario_data.fixed_costs == data.fixed_costs
        @test scenario_data.costs == data.costs
    end

    @testset "StandardSCFLPSubProblem" begin
        cut_strategy = ClassicalCut()
        sp = create_sub_problem(data, cut_strategy)

        # Test struct type and components
        @test sp isa StandardSCFLPSubProblem
        @test length(sp.sub_problems) == data.n_scenarios
        
        # Test each scenario's sub-problem
        for scenario_sp in sp.sub_problems
            @test scenario_sp isa StandardCFLPSubProblem
            @test scenario_sp.model isa Model
            @test length(scenario_sp.fixed_x_constraints) == data.n_facilities

            # Test model constraints for each scenario
            # Test demand satisfaction constraints (=)
            @test num_constraints(scenario_sp.model, AffExpr, MOI.EqualTo{Float64}) == 
                data.n_customers + data.n_facilities  # M demand constraints + N fixed_x constraints

            # Test facility open constraints (<=)
            @test num_constraints(scenario_sp.model, AffExpr, MOI.LessThan{Float64}) == 
                data.n_facilities * data.n_customers + data.n_facilities  # N*M facility open constraints + N capacity constraints

            # Test objective function type
            @test objective_function(scenario_sp.model) isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        end
    end

    @testset "KnapsackSCFLPSubProblem" begin
        cut_strategy = KnapsackCut()
        sp = create_sub_problem(data, cut_strategy)

        # Test struct type and components
        @test sp isa KnapsackSCFLPSubProblem
        @test length(sp.sub_problems) == data.n_scenarios
        
        # Test each scenario's sub-problem
        for scenario_sp in sp.sub_problems
            @test scenario_sp isa KnapsackCFLPSubProblem
            @test scenario_sp.model isa Model
            @test length(scenario_sp.fixed_x_constraints) == data.n_facilities

            # Test model constraints for each scenario
            # Test demand satisfaction constraints (=)
            @test num_constraints(scenario_sp.model, AffExpr, MOI.EqualTo{Float64}) == 
                data.n_customers + data.n_facilities  # M demand constraints + N fixed_x constraints

            # Test knapsack constraints (<=)
            @test num_constraints(scenario_sp.model, AffExpr, MOI.LessThan{Float64}) == 
                data.n_facilities * data.n_customers + data.n_facilities  # N*M facility open constraints + N capacity constraints
            # Test objective function type
            @test objective_function(scenario_sp.model) isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        end
    end
end