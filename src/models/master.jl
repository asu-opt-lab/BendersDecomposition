export Master

# abstract type Master <: AbstractMasterProblem end

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
mutable struct Master <: AbstractMaster
    model::Model
    # obj_value::Float64
    # x_value::Vector{Float64}
    # t_value::Vector{Float64}

    function Master(data::Data)

        @debug "Building Master Problem for CFLP"
    
        model = Model()

        @variable(model, x[1:data.dim_x], Bin)
        @variable(model, t[1:data.dim_t] >= -1e6)

        @objective(model, Min, data.c_x'* x + data.c_t'* t)

        # new(model, 0.0, zeros(data.dim_x), zeros(data.dim_t))
        new(model)#, 0.0, zeros(data.dim_x), zeros(data.dim_t))
    end
end
