export assign_attributes!

const GRB_ENV = Ref{Gurobi.Env}()
function __init__()
    GRB_ENV[] = Gurobi.Env()
    return
end

function assign_attributes!(model::Model, config::Dict{String,Any})
    # Set solver based on config
    # @info config 
    if config["solver"] == "Gurobi"
        set_optimizer(model, () -> Gurobi.Optimizer(GRB_ENV[]))
    elseif config["solver"] == "CPLEX"
        set_optimizer(model, CPLEX.Optimizer)
    else
        error("Unsupported solver: $(config["solver"])")
    end

    # Set solver attributes from config
    for (param, value) in config
        if param != "solver"
            set_optimizer_attribute(model, param, value)
        end
    end

    set_optimizer_attribute(model, MOI.Silent(), true)
end

