module BendersBaseGurobiExt

using BendersBase
using Gurobi
using JuMP
using MathOptInterface
const MOI = MathOptInterface

# Maintain a global Gurobi environment to reuse across models
const GRB_ENV = Ref{Gurobi.Env}()

function _create_gurobi_optimizer(;
    threads::Int = 1,
    mipgap::Float64 = 1e-6,
    feasibilitytol::Float64 = 1e-9,
    intfeastol::Float64 = 1e-9,
    silent::Bool = true,
    kwargs...
)
    function create_optimizer()
        opt = Gurobi.Optimizer(GRB_ENV[])
        return opt
    end
    
    attrs = Pair{Any,Any}[
        "Threads" => threads,
        "MIPGap" => mipgap,
        "FeasibilityTol" => feasibilitytol,
        "IntFeasTol" => intfeastol,
        MOI.Silent() => silent
    ]
    
    for (k, v) in kwargs
        push!(attrs, string(k) => v)
    end
    
    return optimizer_with_attributes(create_optimizer, attrs...)
end

function _create_gurobi_lp_optimizer(;
    threads::Int = 1,
    feasibilitytol::Float64 = 1e-9,
    optimalitytol::Float64 = 1e-9,
    silent::Bool = true,
    kwargs...
)
    function create_optimizer()
        opt = Gurobi.Optimizer(GRB_ENV[])
        return opt
    end
    
    attrs = Pair{Any,Any}[
        "Threads" => threads,
        "FeasibilityTol" => feasibilitytol,
        "OptimalityTol" => optimalitytol,
        MOI.Silent() => silent
    ]
    
    for (k, v) in kwargs
        push!(attrs, string(k) => v)
    end
    
    return optimizer_with_attributes(create_optimizer, attrs...)
end

function __init__()
    # Initialize Gurobi environment
    GRB_ENV[] = Gurobi.Env()
    
    # Register Gurobi as available
    BendersBase.GUROBI_AVAILABLE[] = true
    BendersBase.GUROBI_OPTIMIZER_FACTORY[] = _create_gurobi_optimizer
    BendersBase.GUROBI_LP_OPTIMIZER_FACTORY[] = _create_gurobi_lp_optimizer
    BendersBase.register_solver!(:Gurobi, () -> Gurobi.Optimizer(GRB_ENV[]))
end

end # module
