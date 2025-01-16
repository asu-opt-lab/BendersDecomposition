using Test
using JuMP
using BendersDecomposition

# Helper function to find a constraint that matches a condition
function find_constraint(condition::Function, constraints::Vector{ConstraintRef})
    for con in constraints
        if condition(con)
            return con
        end
    end
    return nothing
end

@testset "SNIP Sub Problem" begin
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
    ψ = [rand(num_nodes) for _ in 1:num_scenarios]

    data = SNIPData(
        num_nodes,
        num_scenarios,
        scenarios,
        D,
        A_minus_D,
        ψ,
        budget
    )

    @testset "StandardSNIPSubProblem" begin
        cut_strategy = ClassicalCut()
        sp = create_sub_problem(data, cut_strategy)

        # Test struct type and components
        @test sp isa StandardSNIPSubProblem
        @test length(sp.sub_problems) == data.num_scenarios

        for k in 1:data.num_scenarios
            sub_k = sp.sub_problems[k]
            @test sub_k isa BendersDecomposition._StandardSNIPSubProblem
            @test sub_k.model isa Model
            @test length(sub_k.fixed_x_constraints) == length(data.D)

            # Get variables
            y = sub_k.model[:y]
            x = sub_k.model[:x]

            # Test variable dimensions
            @test length(y) == data.num_nodes
            @test length(x) == length(data.D)

            # Test y variables are non-negative
            @test all(has_lower_bound(y[i]) for i in 1:data.num_nodes)
            @test all(lower_bound(y[i]) == 0 for i in 1:data.num_nodes)

            # Test objective function
            obj = objective_function(sub_k.model)
            source_node = data.scenarios[k][1]
            @test coefficient(obj, y[source_node]) == 1.0

            # Test destination node constraint
            dest_node = data.scenarios[k][2]
            dest_cons = find_constraint(sub_k.other_constraints) do con
                normalized_coefficient(con, y[dest_node]) == 1 &&
                normalized_rhs(con) == 1
            end
            @test dest_cons !== nothing

            # Test probability propagation constraints
            # For arcs with potential sensors
            for (idx, (from, to, r, q)) in enumerate(data.D)
                # First constraint: y[from] - q * y[to] >= 0
                con1 = find_constraint(sub_k.other_constraints) do con
                    normalized_coefficient(con, y[from]) == 1 &&
                    normalized_coefficient(con, y[to]) == -q &&
                    normalized_rhs(con) == 0
                end
                @test con1 !== nothing

                # Second constraint: y[from] - r * y[to] >= -(r - q) * ψ[k][to] * x[idx]
                con2 = find_constraint(sub_k.other_constraints) do con
                    normalized_coefficient(con, y[from]) == 1 &&
                    normalized_coefficient(con, y[to]) == -r &&
                    normalized_coefficient(con, x[idx]) == (r - q) * data.ψ[k][to] &&
                    normalized_rhs(con) == 0
                end
                @test con2 !== nothing
            end

            # For arcs without sensors
            for (from, to, r) in data.A_minus_D
                con = find_constraint(sub_k.other_constraints) do con
                    normalized_coefficient(con, y[from]) == 1 &&
                    normalized_coefficient(con, y[to]) == -r &&
                    normalized_rhs(con) == 0
                end
                @test con !== nothing
            end
        end
    end
end
