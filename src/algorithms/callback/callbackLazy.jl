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
    lazy_callback(cb_data, master::Master, log::BendersBnBLog, callback::LazyCallback)

Callback function for adding lazy constraints in the branch-and-bound process.
Generates and adds Benders cuts when integer solutions are found.

# Arguments
- `cb_data`: Callback data from the solver
- `master::Master`: The master problem object
- `log::BendersBnBLog`: Log object to record statistics
- `param::BendersBnBParam`: Parameters for the branch-and-bound process
- `callback::LazyCallback`: Configuration for the lazy callback
"""
function lazy_callback(cb_data, master::Master, log::BendersBnBLog, param::BendersBnBParam, callback::LazyCallback)
    status = JuMP.callback_node_status(cb_data, master.model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER

        state = BendersBnBState()
        if solver_name(master.model) == "CPLEX"
            n_count = Ref{CPXINT}()
            ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
            state.node = n_count[]
        end

        state.values[:x] = JuMP.callback_value.(cb_data, master.x)
        state.values[:t] = JuMP.callback_value.(cb_data, master.t)


        state.oracle_time = @elapsed begin
            state.is_in_L, hyperplanes, state.f_x = generate_cuts(callback.oracle, state.values[:x], state.values[:t]; time_limit = param.time_limit)
            cuts = !state.is_in_L ? hyperplanes_to_expression(master.model, hyperplanes, master.x, master.t) : []
            state.num_cuts += length(hyperplanes)
        end

        # Add cuts
        for cut in cuts
            cut_constraint = @build_constraint(0 >= cut)
            MOI.submit(master.model, MOI.LazyConstraint(cb_data), cut_constraint)
        end
        record_node!(log, state, true)

    end
end