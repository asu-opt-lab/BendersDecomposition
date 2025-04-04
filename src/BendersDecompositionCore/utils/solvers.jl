export assign_attributes!

# Single Gurobi environment for all models
const GRB_ENV = Ref{Gurobi.Env}()
function __init__()
    GRB_ENV[] = Gurobi.Env()
end

"""
    assign_attributes!(model::Model, config::Dict{String,Any})

Configure a JuMP model with solver and attributes specified in config dictionary.

# Arguments
- `model::Model`: The JuMP model to configure
- `config::Dict{String,Any}`: Configuration dictionary with solver and parameters
"""
function assign_attributes!(model::Model, config::Dict{String,Any})
    # Set solver based on config
    solver = get(config, "solver", "")
    
    # Configure optimizer
    if solver == "Gurobi"
        set_optimizer(model, () -> Gurobi.Optimizer(GRB_ENV[]))
    elseif solver == "CPLEX" 
        set_optimizer(model, CPLEX.Optimizer)
    else
        error("Unsupported solver: $solver")
    end

    # Apply all non-solver attributes from config
    for (param, value) in config
        param == "solver" && continue
        set_optimizer_attribute(model, param, value)
    end

    # Always set silent mode
    set_optimizer_attribute(model, MOI.Silent(), true)
end
