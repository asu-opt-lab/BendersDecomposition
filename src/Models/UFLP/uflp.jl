

function create_t_variable(model::Model, ::FatKnapsackCut, data::UFLPData)
    M = data.n_customers
    @variable(model, t[1:M] >= -1e6)
end




include("master_problem.jl")
include("sub_problem.jl")
include("dcglp.jl")
include("milp.jl")
