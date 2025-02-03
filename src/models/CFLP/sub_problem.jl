export StandardCFLPSubProblem, KnapsackCFLPSubProblem


abstract type AbstractCFLPSubProblem <: AbstractSubProblem end

"""
    StandardCFLPSubProblem <: AbstractCFLPSubProblem

A mutable struct representing the standard subproblem formulation for the Capacitated Facility Location Problem (CFLP).

# Fields
- `model::Model`: The underlying JuMP optimization model
- `fixed_x_constraints::Vector{ConstraintRef}`: Constraints fixing x as x̂
- `other_constraints::Vector{ConstraintRef}`: Other problem constraints including demand satisfaction and capacity

# Related Functions
    create_sub_problem(data::CFLPData, ::ClassicalCut)
"""
mutable struct StandardCFLPSubProblem <: AbstractCFLPSubProblem
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    other_constraints::Vector{ConstraintRef}
end


struct FacilityKnapsackInfo
    costs::Matrix{Float64}
    demands::Vector{Float64}
    capacity::Vector{Float64}
end

"""
    KnapsackCFLPSubProblem <: AbstractCFLPSubProblem

A mutable struct representing the knapsack-based subproblem formulation for the CFLP.

# Fields
- `model::Model`: The underlying JuMP optimization model
- `fixed_x_constraints::Vector{ConstraintRef}`: Constraints fixing facility opening decisions
- `other_constraints::Vector{ConstraintRef}`: Other problem constraints
- `demand_constraints::Vector{ConstraintRef}`: Customer demand satisfaction constraints for knapsack technique
- `facility_knapsack_info::FacilityKnapsackInfo`: Information for knapsack cut generation

# Related Functions
    create_sub_problem(data::CFLPData, ::KnapsackCut)
"""
mutable struct KnapsackCFLPSubProblem <: AbstractCFLPSubProblem
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    other_constraints::Vector{ConstraintRef}
    demand_constraints::Vector{ConstraintRef}
    facility_knapsack_info::FacilityKnapsackInfo
end



# Specialized create_sub_problem functions
function create_sub_problem(data::CFLPData, ::ClassicalCut)
    @debug "Building Subproblem for CFLP (Standard)"
    model, fixed_x, other, _ = _create_base_sub_problem(data)
    return StandardCFLPSubProblem(model, fixed_x, other)
end

function create_sub_problem(data::CFLPData, ::KnapsackCut)
    @debug "Building Subproblem for CFLP (Knapsack)"
    model, fixed_x, other, demand = _create_base_sub_problem(data)
    facility_knapsack_info = FacilityKnapsackInfo(data.costs, data.demands, data.capacities)
    return KnapsackCFLPSubProblem(model, fixed_x, other, demand, facility_knapsack_info)
end

# Create a base function for common setup logic
function _create_base_sub_problem(data::CFLPData)
    model = Model()
    N, M = data.n_facilities, data.n_customers
    
    # Define variables
    @variable(model, y[1:N, 1:M] >= 0)
    @variable(model, x[1:N])

    # Set objective
    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))

    # Add common constraints
    demand_constraints = @constraint(model, demand_satisfaction[j in 1:M], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:N], sum(data.demands[:] .* y[i,:]) <= data.capacities[i] * x[i])

    # Add initial x constraints
    fixed_x_constraints = @constraint(model, x .== 0)

    # Store other constraints
    other_constraints = Vector{ConstraintRef}()
    append!(other_constraints, model[:demand_satisfaction])
    append!(other_constraints, vec(model[:facility_open]))
    append!(other_constraints, model[:capacity])

    return model, fixed_x_constraints, other_constraints, demand_constraints
end