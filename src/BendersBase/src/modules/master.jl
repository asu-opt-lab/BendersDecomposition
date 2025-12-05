export Master

"""
    Master <: AbstractMaster

A mutable struct representing the master problem in Benders decomposition.

# Fields
- `model::Model`: The underlying JuMP optimization model
- `x::Vector{VariableRef}`: Binary first-stage decision variables
- `t::Vector{VariableRef}`: Auxiliary variables for epigraph constraints
"""
mutable struct Master <: AbstractMaster
    model::Model
    x_tuple::NamedTuple
    x::Vector{VariableRef}
    t::Vector{VariableRef}

    function Master(data::Data; solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9))

        @debug "Building Master Problem for CFLP"
    
        model = Model()

        @variable(model, x[1:data.dim_x], Bin)
        @variable(model, t[1:data.dim_t] >= -1e6)

        @objective(model, Min, data.c_x'* x + data.c_t'* t)

        assign_attributes!(model, solver_param)
        new(model, x, t)
    end

    function Master(problem::AbstractData; customize=customize_master_model!)

        @debug "Building Master module"

        model = Model()

        x_tuple, t = customize(model, problem)
        t = t isa VariableRef ? [t] : t
        
        x = collect(values(x_tuple...))

        new(model, x_tuple, x, t)
    end
end

function customize_master_model!(model::Model, problem::AbstractData)
    throw(UndefError("update customize_master_model! for $(typeof(problem))"))
end
