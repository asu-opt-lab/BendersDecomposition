export LazyCallback

"""
    LazyCallback <: AbstractLazyCallback

Configuration for lazy constraint callbacks in the branch-and-bound process.
Used to dynamically add Benders cuts when integer solutions are found.

# Fields
- `params::EmptyCallbackParam`: Empty parameters for the callback (not used)
- `oracle::AbstractTypicalOracle`: Oracle used to generate Benders cuts (better to be AbstractTypicalOracle as disjunctive oracle at integral node may yield incorrect results.)
"""
struct LazyCallback <: AbstractLazyCallback
    params::EmptyCallbackParam
    oracle::AbstractTypicalOracle
    
    function LazyCallback(oracle::AbstractTypicalOracle)
        new(EmptyCallbackParam(), oracle)
    end
    
    function LazyCallback(data)
        new(EmptyCallbackParam(), ClassicalOracle(data))
    end
end

"""
    lazy_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::LazyCallback)

Callback function for adding lazy constraints in the branch-and-bound process.
Generates and adds Benders cuts when integer solutions are found.

# Arguments
- `cb_data`: Callback data from the solver
- `master_model::Model`: The JuMP master problem model
- `log::BendersBnBLog`: Log object to record statistics
- `param::BendersBnBParam`: Parameters for the branch-and-bound process
- `callback::LazyCallback`: Configuration for the lazy callback
"""
function lazy_callback(cb_data, master_model::Model, log::BendersBnBLog, param::BendersBnBParam, callback::LazyCallback)
    status = JuMP.callback_node_status(cb_data, master_model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER

        state = BendersBnBState()
        if solver_name(master_model) == "CPLEX"
            n_count = Ref{CPXINT}()
            ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
            state.node = n_count[]
        end

        state.values[:x] = JuMP.callback_value.(cb_data, master_model[:x])
        state.values[:t] = JuMP.callback_value.(cb_data, master_model[:t])


        state.oracle_time = @elapsed begin
            state.is_in_L, hyperplanes, state.f_x = generate_cuts(callback.oracle, state.values[:x], state.values[:t]; time_limit = max(get_sec_remaining(log, param), 15))
            cuts = !state.is_in_L ? hyperplanes_to_expression(master_model, hyperplanes, master_model[:x], master_model[:t]) : []
        end
        
        if !isempty(cuts)
            for cut in cuts
                cut_constraint = @build_constraint(0 >= cut)
                MOI.submit(master_model, MOI.LazyConstraint(cb_data), cut_constraint)
                state.num_cuts += 1
            end
        end
        record_node!(log, state, true)

    end
end