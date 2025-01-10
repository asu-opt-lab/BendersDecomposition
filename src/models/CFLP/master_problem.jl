export CFLPMasterProblem


abstract type AbstractCFLPMasterProblem <: AbstractMasterProblem end

mutable struct CFLPMasterProblem <: AbstractCFLPMasterProblem
    model::Model
    var::Dict
    obj_value::Float64
    x_value::Vector{Float64}
    t_value::Float64
end


function create_master_problem(data::CFLPData, cut_strategy::Union{ClassicalCut, KnapsackCut})

    model = Model()

    N = data.n_facilities
    @variable(model, x[1:N], Bin)
    
    t = create_t_variable(model, cut_strategy, data)

    @objective(model, Min, sum(data.fixed_costs .* x) + sum(t))
    @constraint(model, sum(data.capacities[i] * x[i] for i in 1:N) >= sum(data.demands))
    
    return CFLPMasterProblem(model, Dict(:x => x, :t => t), 0.0, zeros(N), 0.0)
end



