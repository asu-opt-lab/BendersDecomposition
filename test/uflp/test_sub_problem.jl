using Test
using JuMP
using Gurobi
using BendersDecomposition

@testset "UFLP Sub Problem" begin
    # Create test data with fixed random seed for reproducibility
    data = UFLPData(
        5, 10,  # 5 facilities, 10 customers
        rand(10) .* 50,     # Customer demands between 0-50
        rand(5) .* 100,     # Facility fixed costs between 0-100 
        rand(5, 10) .* 10,  # Transportation costs between 0-10
    )

    @testset "StandardUFLPSubProblem" begin
        cut_strategy = ClassicalCut()
        sp = create_sub_problem(data, cut_strategy)

        # Test struct type and components
        @test sp isa StandardUFLPSubProblem
        @test sp.model isa Model
        @test length(sp.fixed_x_constraints) == data.n_facilities

        # Test model constraints
        # Should have M demand satisfaction constraints (=) and N facility constraints (=)
        @test num_constraints(sp.model, AffExpr, MOI.EqualTo{Float64}) == 
            data.n_customers + data.n_facilities

        # Should have N*M facility capacity constraints (<=) 
        @test num_constraints(sp.model, AffExpr, MOI.LessThan{Float64}) == 
            data.n_facilities * data.n_customers

        # Test objective function type
        @test objective_function(sp.model) isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
    end

    @testset "KnapsackUFLPSubProblem" begin
        for cut_type in [FatKnapsackCut(), SlimKnapsackCut()]
            sp = create_sub_problem(data, cut_type)

            # Test struct type and components
            @test sp isa KnapsackUFLPSubProblem
            
            # Test cost-demand calculations
            expected_costs = [data.costs[:,j] .* data.demands[j] for j in 1:data.n_customers]
            @test sp.sorted_cost_demands == [sort(costs) for costs in expected_costs]
            
            # Test initialization of data structures
            @test all(isempty(k) for k in values(sp.selected_k))
            @test length(sp.sorted_indices) == data.n_customers
            @test all(length(indices) == data.n_facilities for indices in sp.sorted_indices)
        end
    end
end