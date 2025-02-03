export CFLPMasterProblem

abstract type AbstractCFLPMasterProblem <: AbstractMasterProblem end

"""
    CFLPMasterProblem <: AbstractCFLPMasterProblem

A mutable struct representing the master problem for the Capacitated Facility Location Problem (CFLP).

# Fields
- `model::Model`: The underlying JuMP optimization model
- `var::Dict`: Dictionary storing the problem variables (x, t)
- `obj_value::Float64`: Current objective value of the master problem
- `x_value::Vector{Float64}`: Current values of the integer variables x
- `t_value::Float64`: Current value of the variable t

# Related Functions
    create_master_problem(data::CFLPData, cut_strategy::Union{ClassicalCut, KnapsackCut})
"""
mutable struct CFLPMasterProblem <: AbstractCFLPMasterProblem
    model::Model
    var::Dict
    obj_value::Float64
    x_value::Vector{Float64}
    t_value::Float64
end


function create_master_problem(data::CFLPData, cut_strategy::Union{ClassicalCut, KnapsackCut})
    @debug "Building Master Problem for CFLP"
    
    model = Model()

    N = data.n_facilities
    @variable(model, x[1:N], Bin)
    
    t = create_t_variable(model, cut_strategy, data)

    @objective(model, Min, sum(data.fixed_costs .* x) + sum(t))
    @constraint(model, sum(data.capacities[i] * x[i] for i in 1:N) >= sum(data.demands))
    
    return CFLPMasterProblem(model, Dict(:x => x, :t => t), 0.0, zeros(N), 0.0)
end



