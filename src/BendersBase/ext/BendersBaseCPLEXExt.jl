module BendersBaseCPLEXExt

using BendersBase
using CPLEX
using JuMP
using MathOptInterface
const MOI = MathOptInterface

# Re-export CPLEX types and constants for callback usage
const CPXINT = CPLEX.CPXINT
const CPXcallbackgetinfoint = CPLEX.CPXcallbackgetinfoint
const CPXCALLBACKINFO_NODECOUNT = CPLEX.CPXCALLBACKINFO_NODECOUNT
const CPXCALLBACKINFO_NODEDEPTH = CPLEX.CPXCALLBACKINFO_NODEDEPTH

function _create_cplex_optimizer(; 
    threads::Int = 1,
    epint::Float64 = 1e-9,
    eprhs::Float64 = 1e-9,
    epgap::Float64 = 1e-6,
    silent::Bool = true,
    kwargs...
)
    attrs = Pair{Any,Any}[
        "CPXPARAM_Threads" => threads,
        "CPX_PARAM_EPINT" => epint,
        "CPX_PARAM_EPRHS" => eprhs,
        "CPX_PARAM_EPGAP" => epgap,
        MOI.Silent() => silent
    ]
    
    for (k, v) in kwargs
        push!(attrs, string(k) => v)
    end
    
    return optimizer_with_attributes(CPLEX.Optimizer, attrs...)
end

function _create_cplex_lp_optimizer(;
    threads::Int = 1,
    eprhs::Float64 = 1e-9,
    epopt::Float64 = 1e-9,
    numerical_emphasis::Int = 1,
    silent::Bool = true,
    kwargs...
)
    attrs = Pair{Any,Any}[
        "CPXPARAM_Threads" => threads,
        "CPX_PARAM_EPRHS" => eprhs,
        "CPX_PARAM_EPOPT" => epopt,
        "CPX_PARAM_NUMERICALEMPHASIS" => numerical_emphasis,
        MOI.Silent() => silent
    ]
    
    for (k, v) in kwargs
        push!(attrs, string(k) => v)
    end
    
    return optimizer_with_attributes(CPLEX.Optimizer, attrs...)
end

function __init__()
    # Register CPLEX as available
    BendersBase.CPLEX_AVAILABLE[] = true
    BendersBase.CPLEX_OPTIMIZER_FACTORY[] = _create_cplex_optimizer
    BendersBase.CPLEX_LP_OPTIMIZER_FACTORY[] = _create_cplex_lp_optimizer
    BendersBase.register_solver!(:CPLEX, CPLEX.Optimizer)
end

end # module
