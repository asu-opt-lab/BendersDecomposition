using Test
using JuMP
using BendersDecomposition

@testset "MCNDP Sub Problem" begin
    # Create test data
    num_nodes = 3
    num_arcs = 5
    num_commodities = 2
    
    arcs = [(1,2), (1,3), (2,1), (2,3), (3,1)]
    commodity_paths = [(1,3,10.0), (2,1,15.0)]
    
    data = MCNDPData(
        num_nodes,
        num_arcs,
        num_commodities,
        arcs,
        rand(num_arcs),     # fixed_costs
        rand(num_arcs),     # variable_costs
        rand(num_arcs),     # capacities
        commodity_paths
    )

    @testset "StandardMCNDPSubProblem" begin
        cut_strategy = ClassicalCut()
        sp = create_sub_problem(data, cut_strategy)

        # Test struct type and components
        @test sp isa StandardMCNDPSubProblem
        @test sp.model isa Model
        @test length(sp.fixed_x_constraints) == data.num_arcs
        @test !isempty(sp.other_constraints)

        # Get variables
        x = sp.model[:x]
        y = sp.model[:y]

        # Test variable dimensions
        @test size(y) == (data.num_commodities, data.num_arcs)
        @test length(x) == data.num_arcs

        # Test y variables are non-negative
        @test all(has_lower_bound(y[c,a]) && lower_bound(y[c,a]) == 0 
                 for c in 1:data.num_commodities, a in 1:data.num_arcs)

        # Test objective function
        obj = objective_function(sp.model)
        @test obj isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        for c in 1:data.num_commodities, a in 1:data.num_arcs
            @test coefficient(obj, y[c,a]) ≈ data.variable_costs[a] * data.demands[c][3]
        end

        # Test capacity constraints
        capacity_cons = sp.model[:capacity]
        for a in 1:data.num_arcs
            con = capacity_cons[a]
            @test normalized_coefficient(con, x[a]) ≈ -data.capacities[a]
            for c in 1:data.num_commodities
                @test normalized_coefficient(con, y[c,a]) ≈ data.demands[c][3]
            end
        end

        # Test arc_open constraints
        arc_open_cons = sp.model[:arc_open]
        for c in 1:data.num_commodities, a in 1:data.num_arcs
            con = arc_open_cons[c,a]
            @test normalized_coefficient(con, y[c,a]) == 1
            @test normalized_coefficient(con, x[a]) == -1
        end
    end

    @testset "KnapsackMCNDPSubProblem" begin
        cut_strategy = KnapsackCut()
        sp = create_sub_problem(data, cut_strategy)

        # Test struct type and components
        @test sp isa KnapsackMCNDPSubProblem
        @test sp.model isa Model
        @test length(sp.fixed_x_constraints) == data.num_arcs
        @test !isempty(sp.other_constraints)
        @test size(sp.demand_constraints) == (data.num_commodities, data.num_nodes)
        @test size(sp.b_iv) == (data.num_commodities, data.num_nodes)

        # Test flow conservation constraints
        for c in 1:data.num_commodities
            origin, destination, _ = data.demands[c]
            for i in 1:data.num_nodes
                con = sp.demand_constraints[c,i]
                if i == origin
                    @test sp.b_iv[c,i] == 1
                elseif i == destination
                    @test sp.b_iv[c,i] == -1
                else
                    @test sp.b_iv[c,i] == 0
                end
            end
        end

        # Test that model structure matches StandardMCNDPSubProblem
        x = sp.model[:x]
        y = sp.model[:y]
        @test size(y) == (data.num_commodities, data.num_arcs)
        @test length(x) == data.num_arcs
        @test all(has_lower_bound(y[c,a]) && lower_bound(y[c,a]) == 0 
                 for c in 1:data.num_commodities, a in 1:data.num_arcs)

        # Test objective function
        obj = objective_function(sp.model)
        @test obj isa JuMP.GenericAffExpr{Float64,JuMP.VariableRef}
        for c in 1:data.num_commodities, a in 1:data.num_arcs
            @test coefficient(obj, y[c,a]) ≈ data.variable_costs[a] * data.demands[c][3]
        end
    end
end
