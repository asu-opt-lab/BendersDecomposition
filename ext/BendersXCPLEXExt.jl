module BendersXCPLEXExt

using BendersX
using CPLEX
using JuMP

function _cplex_optimizer(args::Pair...)
    return JuMP.optimizer_with_attributes(
        () -> CPLEX.Optimizer(),
        BendersX._normalize_optimizer_attributes(args...)...,
    )
end

function __init__()
    BendersX._register_cplex_optimizer!(_cplex_optimizer)
    return nothing
end

function BendersX.callback_node_count(cb_data, model::JuMP.Model, ::Val{:CPLEX})
    n_count = Ref{CPXINT}()
    CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
    return Int(n_count[])
end

function BendersX.callback_node_depth(cb_data, model::JuMP.Model, ::Val{:CPLEX})
    node_depth = Ref{CPXINT}()
    CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODEDEPTH, node_depth)
    return Int(node_depth[])
end

end
