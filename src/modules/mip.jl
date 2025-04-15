export Mip

mutable struct Mip <: AbstractMip
    model::Model

    function Mip(data::Data, solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9))

        @debug "Building Mip Problem"
    
        model = Model()

        @variable(model, x[1:data.dim_x], Bin)
        
        assign_attributes!(model, solver_param)

        new(model)
    end
end
