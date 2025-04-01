"""
Helper functions for model creation and manipulation
"""

function create_master_problem(data::AbstractData, cut_strategy::DisjunctiveCut)
    master = create_master_problem(data, cut_strategy.base_cut_strategy)
    relax_integrality(master.model)
    return master
end

function create_t_variable(model::Model, ::CutStrategy, data::AbstractData) 
    @variable(model, t >= -1e6)
    # @variable(model, t >= 0)
end 


function create_sub_problem(data::AbstractData, cut_strategy::DisjunctiveCut)
    sub = create_sub_problem(data, cut_strategy.base_cut_strategy)
    return sub
end
