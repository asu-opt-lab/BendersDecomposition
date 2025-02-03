using Test
using JuMP
using BendersDecomposition

@testset "MCNDP Master Problem" begin
    # Create test data
    num_nodes = 3
    num_arcs = 5
    num_commodities = 2
    
    # Create arc list as tuples of (from_node, to_node)
    arcs = [(1,2), (1,3), (2,1), (2,3), (3,1)]
    
    # Create commodity paths as tuples of (origin, destination, demand)
    commodity_paths = [(1,3,10.0), (2,1,15.0)]
    
    data = MCNDPData(
        num_nodes,           # nodes
        num_arcs,           # arcs
        num_commodities,    # commodities
        arcs,               # list of arcs (from_node, to_node)
        rand(num_arcs),     # fixed_costs
        rand(num_arcs),     # variable_costs
        rand(num_arcs),     # capacities
        commodity_paths     # list of commodities (origin, destination, demand)
    )

    @testset "Classical Cut Strategy" begin
        cut_strategy = ClassicalCut()
        mp = create_master_problem(data, cut_strategy)

        # Test struct type and components
        @test mp isa MCNDPMasterProblem
        @test mp.model isa Model
        @test haskey(mp.var, :x)
        @test haskey(mp.var, :t)

        # Test initial values
        @test mp.obj_value == 0.0
        @test mp.x_value == zeros(data.num_arcs)
        @test mp.t_value == 0.0

        # Test variable dimensions and types
        x = mp.var[:x]
        t = mp.var[:t]
        @test length(x) == data.num_arcs
        @test all(is_binary(x[i]) for i in 1:data.num_arcs)
        @test t isa VariableRef
        @test has_lower_bound(t)
        @test lower_bound(t) == -1e6

        # Test objective function
        obj = objective_function(mp.model)
        @test obj isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        for a in 1:data.num_arcs
            @test coefficient(obj, x[a]) ≈ data.fixed_costs[a]
        end
        @test coefficient(obj, t) == 1.0
    end

    @testset "Knapsack Cut Strategy" begin
        cut_strategy = KnapsackCut()
        mp = create_master_problem(data, cut_strategy)

        # Test struct type and components
        @test mp isa MCNDPMasterProblem
        @test mp.model isa Model
        @test haskey(mp.var, :x)
        @test haskey(mp.var, :t)

        # Test initial values
        @test mp.obj_value == 0.0
        @test mp.x_value == zeros(data.num_arcs)
        @test mp.t_value == 0.0

        # Test variable dimensions and types
        x = mp.var[:x]
        t = mp.var[:t]
        @test length(x) == data.num_arcs
        @test all(is_binary(x[i]) for i in 1:data.num_arcs)
        @test t isa VariableRef
        @test has_lower_bound(t)
        @test lower_bound(t) == -1e6

        # Test objective function
        obj = objective_function(mp.model)
        @test obj isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        for a in 1:data.num_arcs
            @test coefficient(obj, x[a]) ≈ data.fixed_costs[a]
        end
        @test coefficient(obj, t) == 1.0
    end
end
