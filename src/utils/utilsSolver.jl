using Gurobi
using CPLEX

const GRB_ENV = Ref{Gurobi.Env}()
function __init__()
    GRB_ENV[] = Gurobi.Env()
    return
end

"""
    assign_attributes!(model::Model, config::Dict{String,Any})

Assign a solver and solver attributes to `model` from a configuration
dictionary.

The dictionary must contain a `"solver"` entry with value `"Gurobi"` or
`"CPLEX"`. All remaining entries are forwarded as optimizer attributes. The
helper finishes by enabling silent solver output.
"""
function assign_attributes!(model::Model, config::Dict{String,Any})
    # Set solver based on config
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
