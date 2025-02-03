export UFLPMasterProblem

abstract type AbstractUFLPMasterProblem <: AbstractMasterProblem end

"""
    UFLPMasterProblem <: AbstractUFLPMasterProblem

A mutable struct representing the master problem for the Unconstrained Facility Location Problem (UFLP).

# Fields
- `model::Model`: The underlying JuMP optimization model
- `var::Dict`: Dictionary storing the problem variables (x, t)
- `obj_value::Float64`: Current objective value of the master problem
- `x_value::Vector{Float64}`: Current values of the integer variables x
- `t_value::Union{Vector{Float64}, Float64}`: Current values of the variable t

# Related Functions
    create_master_problem(data::UFLPData, cut_strategy::CutStrategy)
"""
mutable struct UFLPMasterProblem <: AbstractUFLPMasterProblem
    model::Model
    var::Dict
    obj_value::Float64
    x_value::Vector{Float64}
    t_value::Union{Vector{Float64}, Float64}
end


function create_master_problem(data::UFLPData, cut_strategy::Union{ClassicalCut, FatKnapsackCut})

    model = Model()

    N = data.n_facilities
    M = data.n_customers
    @variable(model, x[1:N], Bin)
    
    t = create_t_variable(model, cut_strategy, data)

    @objective(model, Min, sum(data.fixed_costs .* x) + sum(t))
    @constraint(model, sum(x) >= 2)

    t_zeros = cut_strategy == FatKnapsackCut() ? zeros(M) : 0.0
    
    return UFLPMasterProblem(model, Dict(:x => x, :t => t), 0.0, zeros(N), t_zeros)
end


