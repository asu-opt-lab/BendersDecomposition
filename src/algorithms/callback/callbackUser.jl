export UserCallbackParam, UserCallback

"""
    UserCallbackParam <: AbstractCallbackParam

Parameters for user callbacks in the branch-and-bound process.

# Fields
- `frequency::Int = 50`: How often to process nodes (every N fractional nodes)
- `node_count::Int = -1`: Only process nodes after this node count (-1 means process all)
- `depth::Int = -1`: Only process nodes with depth >= this value (-1 means process all depths)
"""
Base.@kwdef struct UserCallbackParam <: AbstractCallbackParam
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
    
    function UserCallback(data::Data; params=UserCallbackParam())
        new(params, ClassicalOracle(data))
    end
end

"""
    user_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::UserCallback)

Callback function for adding user cuts in the branch-and-bound process.
Generates and adds Benders cuts at fractional nodes based on the specified frequency and criteria.

# Arguments
- `cb_data`: Callback data from the solver
- `master_model::Model`: The JuMP master problem model
- `log::BendersBnBLog`: Log object to record statistics
- `param::BendersBnBParam`: Parameters for the branch-and-bound process
- `callback::UserCallback`: Configuration for the user callback with parameters controlling when cuts are generated
"""
function user_callback(cb_data, master_model::Model, log::BendersBnBLog, param::BendersBnBParam, callback::UserCallback)
    status = JuMP.callback_node_status(cb_data, master_model)
    
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        log.num_of_fraction_node += 1
        
        # Check if we should process this node based on frequency
        if log.num_of_fraction_node >= callback.params.frequency
            log.num_of_fraction_node = 0
            
            # Get node information if using CPLEX
            process_node = true
            if solver_name(master_model) == "CPLEX" && (callback.params.node_count != -1 || callback.params.depth != -1)
                n_count = Ref{CPXINT}()
                node_depth = Ref{CPXINT}()
                CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
                CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODEDEPTH, node_depth)
                
                # Check if node meets criteria
                if (callback.params.node_count != -1 && n_count[] > callback.params.node_count) || 
                   (callback.params.depth != -1 && node_depth[] < callback.params.depth)
                    process_node = false
                end
            elseif (callback.params.node_count != -1 || callback.params.depth != -1) && solver_name(master_model) != "CPLEX"
                @warn "node_count and depth parameters are not supported for $(solver_name(master_model)) solver"
            end
            
            if process_node
                # Create state and get current variable values
                state = BendersBnBState()
                state.values[:x] = JuMP.callback_value.(cb_data, master_model[:x])
                state.values[:t] = JuMP.callback_value.(cb_data, master_model[:t])
                
                # Generate cuts
                try 
                    state.oracle_time = @elapsed begin
                        state.is_in_L, hyperplanes, state.f_x = generate_cuts(callback.oracle, state.values[:x], state.values[:t]; time_limit = get_sec_remaining(log, param))
                        cuts = !state.is_in_L ? hyperplanes_to_expression(master_model, hyperplanes, master_model[:x], master_model[:t]) : []
                    end
                    # Add cuts if any were generated
                    if !isempty(cuts)
                        for cut in cuts
                            cut_constraint = @build_constraint(0 >= cut)
                            MOI.submit(master_model, MOI.UserCut(cb_data), cut_constraint)
                            state.num_cuts += 1
                        end
                    end
                    
                    # Record node information
                    record_node!(log, state, false)
                catch e
                    if typeof(e) <: TimeLimitException
                        return
                    end
                end
            end
        end
    end
end