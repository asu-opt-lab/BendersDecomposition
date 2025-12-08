export Master

"""
    Master <: AbstractMaster

A mutable struct representing the master problem in Benders decomposition.

# Fields
- `model::Model`: The underlying JuMP optimization model
- `x::Vector{VariableRef}`: Binary first-stage decision variables
- `t::Vector{VariableRef}`: Auxiliary variables for epigraph constraints
- `dim_x::Int`: dimension of `x`
- `dim_t::Int`: dimension of `t`
- `c_x::Vector{Float64}`: objective coefficient vector of `x`
- `c_x::Vector{Float64}`: objective coefficient vector of `t`
"""
mutable struct Master <: AbstractMaster
    model::Model
    x_tuple::NamedTuple
    x::Vector{VariableRef}
    t::Vector{VariableRef}

    dim_x::Int
    dim_t::Int
    c_x::Vector{Float64}
    c_t::Vector{Float64}

    function Master(problem::AbstractData; customize=customize_master_model!)

        @debug "Building Master module"

        model = Model()

        x_tuple, t = customize(model, problem)
        t = t isa VariableRef ? [t] : t
        x = var_from_tuple(x_tuple)

        dim_x = length(x)
        obj = objective_function(model)
        c_x = [coefficient(obj, x[i]) for i in 1:dim_x]
        dim_t = length(t)
        c_t = [coefficient(obj, t[i]) for i in 1:dim_t]

        new(model, x_tuple, x, t, dim_x, dim_t, c_x, c_t)
    end
end

function customize_master_model!(model::Model, problem::AbstractData)
    throw(UndefError("update customize_master_model! for $(typeof(problem))"))
end
