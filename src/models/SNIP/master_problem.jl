export SNIPMasterProblem

abstract type AbstractSNIPMasterProblem <: AbstractMasterProblem end

mutable struct SNIPMasterProblem <: AbstractSNIPMasterProblem
    model::Model
    var::Dict
    obj_value::Float64
    x_value::Vector{Float64}
    t_value::Float64
end

function create_master_problem(data::SNIPData, cut_strategy::Union{ClassicalCut, KnapsackCut})
    model = Model()
    
    # Variables
    @variable(model, x[1:length(data.D)], Bin)  # Binary sensor installation variables
    t = create_t_variable(model, cut_strategy, data)
    
    # Objective
    @objective(model, Min, t)
    
    # Sensor budget constraint
    @constraint(model, sum(x) <= data.budget)
    
    return SNIPMasterProblem(model, Dict(:x => x, :t => t), 0.0, zeros(length(data.D)), 0.0)
end 