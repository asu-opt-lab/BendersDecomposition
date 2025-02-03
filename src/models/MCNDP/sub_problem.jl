export StandardMCNDPSubProblem, KnapsackMCNDPSubProblem

abstract type AbstractMCNDPSubProblem <: AbstractSubProblem end

"""
    StandardMCNDPSubProblem <: AbstractMCNDPSubProblem

A mutable struct representing the standard subproblem formulation for the Multi-Commodity Network Design Problem (MCNDP).

# Fields
- `model::Model`: The underlying JuMP optimization model
- `fixed_x_constraints::Vector{ConstraintRef}`: Constraints fixing arc selection decisions (x)
- `other_constraints::Vector{ConstraintRef}`: Problem constraints including flow conservation and capacity limits

# Related Functions
    create_sub_problem(data::MCNDPData, ::ClassicalCut)
"""
mutable struct StandardMCNDPSubProblem <: AbstractMCNDPSubProblem
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    other_constraints::Vector{ConstraintRef}
end

"""
    KnapsackMCNDPSubProblem <: AbstractMCNDPSubProblem

A mutable struct representing the knapsack-based subproblem formulation for the MCNDP.

# Fields
- `model::Model`: The underlying JuMP optimization model
- `fixed_x_constraints::Vector{ConstraintRef}`: Constraints fixing arc selection decisions
- `other_constraints::Vector{ConstraintRef}`: General problem constraints
- `demand_constraints::Matrix{ConstraintRef}`: Flow conservation constraints for each commodity-node pair
- `data::MCNDPData`: Problem instance data
- `b_iv::Matrix{Float64}`: Right-hand side values for flow conservation constraints

# Related Functions
    create_sub_problem(data::MCNDPData, ::KnapsackCut)
"""
mutable struct KnapsackMCNDPSubProblem <: AbstractMCNDPSubProblem
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    other_constraints::Vector{ConstraintRef}
    demand_constraints::Matrix{ConstraintRef}
    data::MCNDPData
    b_iv::Matrix{Float64}
end

# Specialized create_sub_problem functions
function create_sub_problem(data::MCNDPData, ::ClassicalCut)
    @debug "Building Subproblem for MCNDP (Standard)"
    model, fixed_x, other, _, _ = _create_base_sub_problem(data)
    return StandardMCNDPSubProblem(model, fixed_x, other)
end

function create_sub_problem(data::MCNDPData, ::KnapsackCut)
    @debug "Building Subproblem for MCNDP (Knapsack Cut)"
    model, fixed_x, other, demand_constraints, b_iv = _create_base_sub_problem(data)
    return KnapsackMCNDPSubProblem(model, fixed_x, other, demand_constraints, data, b_iv)
end

# Create a base function for common setup logic
function _create_base_sub_problem(data::MCNDPData)
    model = Model()
    @variable(model, x[1:data.num_arcs])  # Binary decision variables for arc selection
    @variable(model, y[1:data.num_commodities, 1:data.num_arcs] >= 0)  # Flow variables
    
    # Objective function
    @objective(model, Min, 
        sum(data.variable_costs[a] * data.demands[c][3] * y[c,a] 
            for a in 1:data.num_arcs, c in 1:data.num_commodities)
    )
    
    # Capacity constraints
    @constraint(model, capacity[a in 1:data.num_arcs],
        sum(data.demands[c][3] * y[c,a] for c in 1:data.num_commodities) <= data.capacities[a] * x[a]
    )
    
    @constraint(model, arc_open[c in 1:data.num_commodities, a in 1:data.num_arcs], y[c,a] <= x[a])

    demand_constraints = Matrix{ConstraintRef}(undef, data.num_commodities, data.num_nodes)
    b_iv = Matrix{Float64}(undef, data.num_commodities, data.num_nodes)
    # Flow conservation constraints
    for c in 1:data.num_commodities # for each commodity
        origin, destination, _ = data.demands[c]
        
        for i in 1:data.num_nodes # for each node
            # Calculate inflow and outflow
            inflow = sum(y[c,a] for a in 1:data.num_arcs if data.arcs[a][2] == i)
            outflow = sum(y[c,a] for a in 1:data.num_arcs if data.arcs[a][1] == i)
            
            # Node balance constraint depends on whether it's origin, destination or intermediate
            if i == origin
                demand_constraints[c,i] = @constraint(model, outflow - inflow == 1)
                b_iv[c,i] = normalized_rhs(demand_constraints[c,i])
            elseif i == destination
                demand_constraints[c,i] = @constraint(model, outflow - inflow == -1)
                b_iv[c,i] = normalized_rhs(demand_constraints[c,i])
            else
                demand_constraints[c,i] = @constraint(model, outflow - inflow == 0)
                b_iv[c,i] = normalized_rhs(demand_constraints[c,i])
            end
        end
    end

    # Add initial x constraints
    fixed_x_constraints = @constraint(model, x .== 0)
    other_constraints = Vector{ConstraintRef}()
    append!(other_constraints, model[:arc_open])
    append!(other_constraints, model[:capacity])
    append!(other_constraints, demand_constraints)

    return model, fixed_x_constraints, other_constraints, demand_constraints, b_iv
end