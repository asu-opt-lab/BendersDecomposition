export SCFLPMasterProblem

abstract type AbstractSCFLPMasterProblem <: AbstractMasterProblem end

"""
    SCFLPMasterProblem <: AbstractSCFLPMasterProblem

A mutable struct representing the master problem for the Stochastic Capacitated Facility Location Problem (SCFLP).

# Fields
- `model::Model`: The underlying JuMP optimization model
- `var::Dict`: Dictionary storing the problem variables (x, t)
- `obj_value::Float64`: Current objective value of the master problem
- `x_value::Vector{Float64}`: Current values of the integer variables x
- `t_value::Vector{Float64}`: Current values of the variables t for each scenario

# Related Functions
    create_master_problem(data::SCFLPData, cut_strategy::Union{ClassicalCut, KnapsackCut})
"""
mutable struct SCFLPMasterProblem <: AbstractSCFLPMasterProblem
    model::Model
    var::Dict
    obj_value::Float64
    x_value::Vector{Float64}
    t_value::Vector{Float64}
end

function create_master_problem(data::SCFLPData, cut_strategy::Union{ClassicalCut, KnapsackCut})
    model = Model()

    N = data.n_facilities
    @variable(model, x[1:N], Bin)
    
    t = create_t_variable(model, cut_strategy, data)

    @objective(model, Min, sum(data.fixed_costs .* x) + sum(t)/data.n_scenarios)

    # Add capacity constraint for maximum demand scenario
    max_demand = maximum(sum(demands) for demands in data.demands)
    @constraint(model, sum(data.capacities[i] * x[i] for i in 1:N) >= max_demand)
    
    return SCFLPMasterProblem(model, Dict(:x => x, :t => t), 0.0, zeros(N), zeros(data.n_scenarios))
end

function create_t_variable(model::Model, cut_strategy::Union{ClassicalCut, KnapsackCut}, data::SCFLPData)
    @variable(model, t[1:data.n_scenarios] >= -1e6)
end
