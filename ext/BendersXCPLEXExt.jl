module BendersXCPLEXExt

using BendersX
using CPLEX
using JuMP

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
