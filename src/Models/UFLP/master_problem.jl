export UFLPMasterProblem


abstract type AbstractUFLPMasterProblem <: AbstractMasterProblem end

mutable struct UFLPMasterProblem <: AbstractUFLPMasterProblem
    model::Model
    var::Dict
    obj_value::Float64
    x_value::Vector{Float64}
end


function create_master_problem(data::UFLPData, cut_strategy::CutGenerationStrategy)

    model = Model()

    N = data.n_facilities
    @variable(model, x[1:N], Bin)
    
    t = create_t_variable(model, cut_strategy, data)

    @objective(model, Min, sum(data.fixed_costs .* x) + sum(t))
    @constraint(model, sum(x) >= 2)
    
    return UFLPMasterProblem(model, Dict(:x => x, :t => t), 0.0, zeros(N))
end

# Helper functions for variable creation
function create_t_variable(model::Model, ::CutGenerationStrategy, data::UFLPData) 
    @variable(model, t >= -1e6)
end

function create_t_variable(model::Model, ::FatKnapsackCut, data::UFLPData)
    M = data.n_customers
    @variable(model, t[1:M] >= -1e6)
end
