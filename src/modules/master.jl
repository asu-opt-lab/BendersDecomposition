export Master

mutable struct Master <: AbstractMaster
    model::Model
    
    function Master(data::Data; solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9))

        @debug "Building Master Problem for CFLP"
    
        model = Model()

        @variable(model, x[1:data.dim_x], Bin)
        @variable(model, t[1:data.dim_t] >= -1e6)

        @objective(model, Min, data.c_x'* x + data.c_t'* t)

        assign_attributes!(model, solver_param)
        
        new(model)
    end
end
