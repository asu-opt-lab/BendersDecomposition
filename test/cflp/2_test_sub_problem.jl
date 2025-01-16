using Test
using JuMP
using Gurobi
using BendersDecomposition

@testset "CFLP Sub Problem" begin
    # Create test data with fixed random seed for reproducibility
    data = CFLPData(
        5, 10,  # 5 facilities, 10 customers
        rand(5) .* 200,     # Random capacities between 0-200
        rand(10) .* 50,     # Customer demands between 0-50
        rand(5) .* 100,     # Facility fixed costs between 0-100
        rand(5, 10) .* 10,  # Transportation costs between 0-10
    )

    @testset "StandardCFLPSubProblem" begin
        cut_strategy = ClassicalCut()
        sp = create_sub_problem(data, cut_strategy)

        # Test struct type and components
        @test sp isa StandardCFLPSubProblem
        @test sp.model isa Model
        @test length(sp.fixed_x_constraints) == data.n_facilities
        @test !isempty(sp.other_constraints)

        # Test model constraints
        # Test demand satisfaction constraints (=)
        @test num_constraints(sp.model, AffExpr, MOI.EqualTo{Float64}) == 
            data.n_customers + data.n_facilities  # M demand constraints + N fixed_x constraints

        # Test facility open constraints (<=)
        @test num_constraints(sp.model, AffExpr, MOI.LessThan{Float64}) == 
            data.n_facilities * data.n_customers + data.n_facilities  # N*M facility open constraints + N capacity constraints

        # Test objective function type
        @test objective_function(sp.model) isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
    end

    @testset "KnapsackCFLPSubProblem" begin
        cut_strategy = KnapsackCut()
        sp = create_sub_problem(data, cut_strategy)

        # Test struct type and components
        @test sp isa KnapsackCFLPSubProblem
        @test sp.model isa Model
        @test length(sp.fixed_x_constraints) == data.n_facilities
        @test !isempty(sp.other_constraints)
        @test length(sp.demand_constraints) == data.n_customers

        # Test facility knapsack info
        @test sp.facility_knapsack_info isa BendersDecomposition.FacilityKnapsackInfo
        @test size(sp.facility_knapsack_info.costs) == (data.n_facilities, data.n_customers)
        @test length(sp.facility_knapsack_info.demands) == data.n_customers
        @test length(sp.facility_knapsack_info.capacity) == data.n_facilities

        # Test model constraints (same structure as StandardCFLPSubProblem)
        @test num_constraints(sp.model, AffExpr, MOI.EqualTo{Float64}) == 
            data.n_customers + data.n_facilities
        @test num_constraints(sp.model, AffExpr, MOI.LessThan{Float64}) == 
            data.n_facilities * data.n_customers + data.n_facilities
    end

end