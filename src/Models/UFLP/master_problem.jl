export UFLPMasterProblem


abstract type AbstractUFLPMasterProblem <: AbstractMasterProblem end

mutable struct UFLPMasterProblem <: AbstractUFLPMasterProblem
    model::Model
    var::Dict
    obj_value::Float64
    x_value::Vector{Float64}
    t_value::Vector{Float64}
end


function create_master_problem(data::UFLPData, cut_strategy::Union{ClassicalCut,FatKnapsackCut}
)

    model = Model()

    N = data.n_facilities
    M = data.n_customers
    @variable(model, x[1:N], Bin)
    
    t = create_t_variable(model, cut_strategy, data)

    @objective(model, Min, sum(data.fixed_costs .* x) + sum(t))
    @constraint(model, sum(x) >= 2)
    
    return UFLPMasterProblem(model, Dict(:x => x, :t => t), 0.0, zeros(N), zeros(M))
end


