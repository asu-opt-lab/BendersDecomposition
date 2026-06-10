
"""
    LazyCallback <: AbstractLazyCallback

Configuration for lazy constraint callbacks in the branch-and-bound process.
Used to dynamically add Benders cuts when integer solutions are found.

# Fields
- `oracle::AbstractTypicalOracle`: Oracle used to generate Benders cuts

# Notes
- It is recommended to use `AbstractTypicalOracle` rather than disjunctive oracles, as disjunctive oracles at integral nodes may yield incorrect results due to the nature of disjunctive programming.
"""
struct LazyCallback <: AbstractLazyCallback
    oracle::AbstractTypicalOracle
    
    function LazyCallback(oracle::AbstractTypicalOracle)
        new(oracle)
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

# Notes
- Solver-specific callback metadata such as node counts is resolved through
  [`callback_node_count`](@ref).
- CPLEX support for this accessor is provided by `BendersXCPLEXExt` when
  `CPLEX.jl` is loaded.
"""
function record_node_count!(state::BendersBnBState, node_count::Union{Nothing,Int})
    if !isnothing(node_count)
        state.node = node_count
    end
    return state
end

function lazy_callback(cb_data, master::Master, log::BendersBnBLog, param::BendersBnBParam, callback::LazyCallback)
    status = JuMP.callback_node_status(cb_data, master.model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER

        state = BendersBnBState()
        record_node_count!(state, callback_node_count(cb_data, master.model))

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
