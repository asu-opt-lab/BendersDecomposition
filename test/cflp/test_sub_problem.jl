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

end