export StandardUFLPSubProblem, KnapsackUFLPSubProblem

abstract type AbstractUFLPSubProblem <: AbstractSubProblem end

"""
    StandardUFLPSubProblem <: AbstractUFLPSubProblem

A mutable struct representing the subproblem for the Unconstrained Facility Location Problem (UFLP) with a classical cut.

# Fields
- `model::Model`: The underlying JuMP optimization model
- `fixed_x_constraints::Vector{ConstraintRef}`: Constraints fixing the x variables to 0

# Related Functions
    create_sub_problem(data::UFLPData, cut_strategy::ClassicalCut)
"""
mutable struct StandardUFLPSubProblem <: AbstractUFLPSubProblem
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
end

"""
    KnapsackUFLPSubProblem <: AbstractUFLPSubProblem

A mutable struct representing the subproblem for the Unconstrained Facility Location Problem (UFLP) with a knapsack cut.

# Fields
- `sorted_cost_demands::Vector{Vector{Float64}}`: Sorted cost demands for each customer
- `sorted_indices::Vector{Vector{Int}}`: Sorted indices for each customer
- `selected_k::Dict`: Dictionary storing the selected k values for each customer

# Related Functions
    create_sub_problem(data::UFLPData, cut_strategy::Union{FatKnapsackCut, SlimKnapsackCut})
"""
mutable struct KnapsackUFLPSubProblem <: AbstractUFLPSubProblem
    sorted_cost_demands::Vector{Vector{Float64}}
    sorted_indices::Vector{Vector{Int}}
    selected_k::Dict
end


function create_sub_problem(data::UFLPData, ::ClassicalCut)
    
    model = Model()
    
    @debug "Building Subproblem for UFLP (Standard)"

    # Extract problem dimensions
    N, M = data.n_facilities, data.n_customers
    
    # Define variables
    @variable(model, y[1:N, 1:M] >= 0)
    @variable(model, x[1:N])

    # Set objective
    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))

    # Add constraints
    @constraint(model, demand_satisfaction[j in 1:M], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)

    # Add constraints to fix x variables
    fixed_x_constraints = @constraint(model, x .== 0)  # Initial values, will be updated later

    # Create and return the StandardUFLPSubProblem struct
    return StandardUFLPSubProblem(model, fixed_x_constraints)
end

function create_sub_problem(data::UFLPData, cut_strategy::Union{FatKnapsackCut, SlimKnapsackCut})
    
    @debug "Building Subproblem for UFLP ($(typeof(cut_strategy)))"
    cost_demands = [data.costs[:,j] .* data.demands[j] for j in 1:data.n_customers]
    sorted_indices = [sortperm(cost_demands[j]) for j in 1:data.n_customers]
    sorted_cost_demands = [cost_demands[j][sorted_indices[j]] for j in 1:data.n_customers]
    selected_k = Dict(j => [] for j in 1:data.n_customers)
    return KnapsackUFLPSubProblem(sorted_cost_demands, sorted_indices, selected_k)
end
