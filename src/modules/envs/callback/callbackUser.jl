
"""
    UserCallbackParam <: AbstractUserCallbackParam

Parameters for user callbacks in the branch-and-bound process.

# Fields
- `frequency::Int = 50`: How often to process nodes (every N fractional nodes)
- `node_count::Int = -1`: Minimum solver node count to process (-1 means no node-count threshold)
- `depth::Int = -1`: Minimum solver node depth to process (-1 means no depth threshold)
"""
Base.@kwdef struct UserCallbackParam <: AbstractUserCallbackParam
    frequency::Int = 50
    node_count::Int = -1
    depth::Int = -1
end


"""
    UserCallback <: AbstractUserCallback

Configuration for user cut callbacks in the branch-and-bound process.
Used to dynamically add Benders cuts at fractional nodes.

# Fields
- `params::UserCallbackParam`: Parameters controlling when cuts are generated
- `oracle::AbstractOracle`: Oracle used to generate Benders cuts
"""
struct UserCallback <: AbstractUserCallback
    params::UserCallbackParam
    oracle::AbstractOracle
    
    function UserCallback(oracle::AbstractOracle; params=UserCallbackParam())
        new(params, oracle)
    end
end

"""
    user_callback(cb_data, master::Master, log::BendersBnBLog, callback::UserCallback)

Callback function for adding user cuts in the branch-and-bound process.
Generates and adds Benders cuts at fractional nodes based on the specified frequency and criteria.

# Arguments
- `cb_data`: Callback data from the solver
- `master::Master`: The master module
- `log::BendersBnBLog`: Log object to record statistics
- `param::BendersBnBParam`: Parameters for the branch-and-bound process
- `callback::UserCallback`: Configuration for the user callback with parameters controlling when cuts are generated

# Notes
- Solver-specific callback metadata such as node count and depth is resolved through
  [`callback_node_count`](@ref) and [`callback_node_depth`](@ref).
- CPLEX support for these accessors is provided by `BendersXCPLEXExt` when `CPLEX.jl`
  is loaded.
- If callback metadata is unsupported by the active solver, the fallback metadata accessor emits a warning and the corresponding threshold is not applied.
"""

function user_callback(cb_data, master::Master, log::BendersBnBLog, param::BendersBnBParam, callback::UserCallback)
    status = JuMP.callback_node_status(cb_data, master.model)
    
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        log.num_of_fraction_node += 1
        
        # Check if we should process this node based on frequency
        if log.num_of_fraction_node >= callback.params.frequency
            log.num_of_fraction_node = 0
            
            node_count = callback.params.node_count == -1 ? nothing :
                callback_node_count(cb_data, master.model)
            node_depth = callback.params.depth == -1 ? nothing :
                callback_node_depth(cb_data, master.model)

            process_node =
                (isnothing(node_count) || node_count >= callback.params.node_count) &&
                (isnothing(node_depth) || node_depth >= callback.params.depth)
            
            if process_node
                # Create state and get current variable values
                state = BendersBnBState()
                state.values[:x] = JuMP.callback_value.(cb_data, master.x)
                state.values[:t] = JuMP.callback_value.(cb_data, master.t)
                
                # Generate cuts
                state.oracle_time = @elapsed begin
                    state.is_in_L, hyperplanes, state.f_x = generate_cuts(callback.oracle, state.values[:x], state.values[:t]; time_limit = get_sec_remaining(log, param))
                    cuts = !state.is_in_L ? hyperplanes_to_expression(master.model, hyperplanes, master.x, master.t) : []
                    state.num_cuts += length(hyperplanes)
                end

                # Add cuts
                for cut in cuts
                    cut_constraint = @build_constraint(0 >= cut)
                    MOI.submit(master.model, MOI.UserCut(cb_data), cut_constraint)
                end
                
                # Record node information
                record_node!(log, state, false)
            end
        end
    end
end
