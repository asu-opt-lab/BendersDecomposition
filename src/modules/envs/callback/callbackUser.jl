
"""
    UserCallbackParam <: AbstractUserCallbackParam

Parameters for user callbacks in the branch-and-bound process.

# Fields
- `frequency::Int = 50`: How often to process nodes (every N fractional nodes)
- `node_count::Int = -1`: Minimum node count threshold; only process nodes with `node_count >= threshold` (-1 means process all)
- `depth::Int = -1`: Only process nodes with depth >= this value (-1 means process all depths)
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
- Solver-specific callback metadata such as node count and depth is resolved through package extensions.
"""
function _evaluate_user_callback_filters(
    params::UserCallbackParam;
    node_count::Union{Nothing,Int} = nothing,
    node_depth::Union{Nothing,Int} = nothing,
)
    missing_node_count = params.node_count != -1 && isnothing(node_count)
    missing_depth = params.depth != -1 && isnothing(node_depth)

    process_node = true
    if params.node_count != -1 && !isnothing(node_count) && node_count < params.node_count
        process_node = false
    end
    if params.depth != -1 && !isnothing(node_depth) && node_depth < params.depth
        process_node = false
    end

    return (; process_node, missing_node_count, missing_depth)
end

function _warn_ignored_user_callback_filters(
    solver::AbstractString;
    missing_node_count::Bool = false,
    missing_depth::Bool = false,
)
    ignored_filters = String[]
    missing_node_count && push!(ignored_filters, "node_count")
    missing_depth && push!(ignored_filters, "depth")
    isempty(ignored_filters) && return nothing

    warning_id = if missing_node_count && missing_depth
        :user_callback_missing_node_count_depth
    elseif missing_node_count
        :user_callback_missing_node_count
    else
        :user_callback_missing_depth
    end
    ignored_filter_list = join(ignored_filters, ", ")
    message = "Ignoring unsupported user callback filter(s) $(ignored_filter_list) for $(solver) solver because the required callback metadata is unavailable."
    @warn message maxlog=1 _id=warning_id
    return nothing
end

function user_callback(cb_data, master::Master, log::BendersBnBLog, param::BendersBnBParam, callback::UserCallback)
    status = JuMP.callback_node_status(cb_data, master.model)
    
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        log.num_of_fraction_node += 1
        
        # Check if we should process this node based on frequency
        if log.num_of_fraction_node >= callback.params.frequency
            log.num_of_fraction_node = 0
            
            node_count = callback.params.node_count == -1 ? nothing : callback_node_count(cb_data, master.model)
            node_depth = callback.params.depth == -1 ? nothing : callback_node_depth(cb_data, master.model)
            filter_result = _evaluate_user_callback_filters(callback.params; node_count=node_count, node_depth=node_depth)
            _warn_ignored_user_callback_filters(
                solver_name(master.model);
                missing_node_count=filter_result.missing_node_count,
                missing_depth=filter_result.missing_depth,
            )
            
            if filter_result.process_node
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
