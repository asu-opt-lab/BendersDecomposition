export MCNDPMasterProblem

abstract type AbstractMCNDPMasterProblem <: AbstractMasterProblem end

"""
    MCNDPMasterProblem <: AbstractMCNDPMasterProblem

A mutable struct representing the master problem for the Multi-Commodity Network Design Problem (MCNDP).

# Fields
- `model::Model`: The underlying JuMP optimization model
- `var::Dict`: Dictionary storing the problem variables (x, t)
- `obj_value::Float64`: Current objective value of the master problem
- `x_value::Vector{Float64}`: Current values of the integer variables x for arc selection
- `t_value::Float64`: Current value of the variable t

# Related Functions
    create_master_problem(data::MCNDPData, cut_strategy::Union{ClassicalCut, KnapsackCut})
"""
mutable struct MCNDPMasterProblem <: AbstractMCNDPMasterProblem
    model::Model
    var::Dict
    obj_value::Float64
    x_value::Vector{Float64}
    t_value::Float64
end

function create_master_problem(data::MCNDPData, cut_strategy::Union{ClassicalCut, KnapsackCut})
    model = Model()

    @variable(model, x[1:data.num_arcs], Bin)  # Binary decision variables for arc selection
    t = create_t_variable(model, cut_strategy, data)
    
    # Objective function
    @objective(model, Min, 
        sum(data.fixed_costs[a] * x[a] for a in 1:data.num_arcs) + t
    )

    return MCNDPMasterProblem(model, Dict(:x => x, :t => t), 0.0, zeros(data.num_arcs), 0.0)
end

function create_t_variable(model::Model, cut_strategy::Union{ClassicalCut, KnapsackCut}, data::MCNDPData)
    @variable(model, t >= -1e6)
end
