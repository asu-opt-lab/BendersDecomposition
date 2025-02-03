export SNIPMasterProblem

abstract type AbstractSNIPMasterProblem <: AbstractMasterProblem end

"""
    SNIPMasterProblem <: AbstractSNIPMasterProblem

A mutable struct representing the master problem for the Sensor Network Installation Problem (SNIP).

# Fields
- `model::Model`: The underlying JuMP optimization model
- `var::Dict`: Dictionary storing the problem variables (x, t)
- `obj_value::Float64`: Current objective value of the master problem
- `x_value::Vector{Float64}`: Current values of the integer variables x
- `t_value::Vector{Float64}`: Current values of the variable t

# Related Functions
    create_master_problem(data::SNIPData, cut_strategy::Union{ClassicalCut, KnapsackCut})
"""
mutable struct SNIPMasterProblem <: AbstractSNIPMasterProblem
    model::Model
    var::Dict
    obj_value::Float64
    x_value::Vector{Float64}
    t_value::Vector{Float64}
end

function create_master_problem(data::SNIPData, cut_strategy::ClassicalCut)
    model = Model()
    
    # Variables
    @variable(model, x[1:length(data.D)], Bin)  # Binary sensor installation variables
    t = create_t_variable(model, cut_strategy, data)
    
    # Objective
    @objective(model, Min, sum(data.scenarios[k][3] * t[k] for k in 1:data.num_scenarios))
    
    # Sensor budget constraint
    @constraint(model, sum(x) <= data.budget)
    
    return SNIPMasterProblem(model, Dict(:x => x, :t => t), 0.0, zeros(length(data.D)), zeros(data.num_scenarios))
end 


function create_t_variable(model::Model, cut_strategy::ClassicalCut, data::SNIPData)
    @variable(model, t[1:data.num_scenarios] >= -1e6)
end