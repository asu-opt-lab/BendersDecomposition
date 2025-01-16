using Test
using JuMP
using BendersDecomposition

@testset "SNIP Master Problem" begin
    # Create test data
    num_nodes = 5
    num_scenarios = 3
    budget = 2.0
    
    # Create scenarios as (source, target, probability)
    scenarios = [
        (1, 5, 0.4),
        (2, 4, 0.3),
        (3, 5, 0.3)
    ]
    
    # Create D as (from_node, to_node, r, q)
    D = [
        (1, 2, 1.0, 0.9),
        (1, 3, 2.0, 0.8),
        (2, 3, 1.0, 0.9),
        (2, 4, 2.0, 0.8),
        (3, 4, 1.0, 0.9),
        (3, 5, 2.0, 0.8),
        (4, 5, 1.0, 0.9)
    ]
    
    # Create A_minus_D as (from_node, to_node, r)
    A_minus_D = [
        (1, 4, 3.0),
        (1, 5, 4.0),
        (2, 5, 3.0),
        (4, 1, 3.0),
        (5, 1, 4.0),
        (5, 2, 3.0)
    ]
    
    # Create ψ matrix
    ψ = [rand(length(D)) for _ in 1:num_scenarios]

    data = SNIPData(
        num_nodes,
        num_scenarios,
        scenarios,
        D,
        A_minus_D,
        ψ,
        budget
    )

    @testset "Classical Cut Strategy" begin
        cut_strategy = ClassicalCut()
        mp = create_master_problem(data, cut_strategy)

        # Test struct type and components
        @test mp isa SNIPMasterProblem
        @test mp.model isa Model
        @test haskey(mp.var, :x)
        @test haskey(mp.var, :t)

        # Test initial values
        @test mp.obj_value == 0.0
        @test mp.x_value == zeros(length(data.D))
        @test mp.t_value == zeros(data.num_scenarios)

        # Test variable dimensions and types
        x = mp.var[:x]
        t = mp.var[:t]
        @test length(x) == length(data.D)
        @test length(t) == data.num_scenarios
        @test all(is_binary(x[i]) for i in 1:length(data.D))
        @test all(has_lower_bound(t[k]) for k in 1:data.num_scenarios)
        @test all(lower_bound(t[k]) == -1e6 for k in 1:data.num_scenarios)

        # Test objective function
        obj = objective_function(mp.model)
        @test obj isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        for k in 1:data.num_scenarios
            @test coefficient(obj, t[k]) ≈ data.scenarios[k][3]
        end

        # Test budget constraint
        budget_cons = all_constraints(mp.model, AffExpr, MOI.LessThan{Float64})
        @test length(budget_cons) == 1
        @test normalized_rhs(budget_cons[1]) == data.budget
        @test all(normalized_coefficient(budget_cons[1], x[i]) == 1 
                 for i in 1:length(data.D))
    end
end
