struct GenericCallback <: AbstractCallback 
    time_limit::Union{Float64,Nothing}
    gap_tolerance::Union{Float64}
    lazy_callback::Function
    user_callback::Union{Function, Nothing}
    verbose::Bool
end

function solve!(env::AbstractBendersEnv, solution_strategy::GenericCallback, cut_strategy::AbstractCutStrategy)
    # Create a closure that captures env and cut_strategy
    function callback_wrapper(cb_data)
        try
            @debug "Starting lazy constraint callback"
            lazy_callback(cb_data, env, cut_strategy)
        catch e
            @error "Error in lazy constraint callback" exception=(e, catch_backtrace())
            rethrow(e)
        end
    end
    return _solve!(env, solution_strategy, cut_strategy, callback_wrapper)
end

function lazy_callback(cb_data, env::AbstractBendersEnv, cut_strategy::AbstractCutStrategy)
    status = JuMP.callback_node_status(cb_data, env.master.model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER 
        try
            n_count = Ref{CPXINT}()
            ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
            @debug "Processing node" node_count=n_count[]

            # Get current solution values
            x_values = JuMP.callback_value.(cb_data, env.master.variables[:integer_variables])
            t_values = JuMP.callback_value.(cb_data, env.master.variables[:continuous_variables])
            
            # Solve subproblem and generate cuts
            solve_sub!(env.sub, x_values)
            cuts, sub_obj_value = generate_cuts(env, cut_strategy)
            add_cuts!(env, cuts, sub_obj_value, cb_data)
            
            @debug "Added cuts" num_cuts=length(cuts) sub_obj_value=sub_obj_value
        catch e
            @error "Error in lazy constraint processing" exception=(e, catch_backtrace())
            rethrow(e)
        end
    end
end

function _solve!(env::AbstractBendersEnv, solution_strategy::GenericCallback, cut_strategy::AbstractCutStrategy, lazy_callback)
    start_time = time()

    # Root node preprocessing
    preprocessing_time = @elapsed begin
        df_root_node_preprocessing = root_node_preprocessing!(env, cut_strategy, solution_strategy)
    end
    
    # Adjust time limit after preprocessing
    remaining_time = max(0.0, solution_strategy.time_limit - preprocessing_time)
    
    # Configure solver
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)
    set_time_limit_sec(env.master.model, remaining_time)
    set_optimizer_attribute(env.master.model, MOI.Silent(), !solution_strategy.verbose)
    set_optimizer_attribute(env.master.model, "CPXPARAM_MIP_Display", solution_strategy.verbose ? 3 : 0)

    # Solve the model
    solve_time = @elapsed begin
        JuMP.optimize!(env.master.model)
    end

    # Collect results
    results = Dict(
        :node_count => JuMP.node_count(env.master.model),
        :elapsed_time => time() - start_time,
        :objective_bound => JuMP.objective_bound(env.master.model),
        :objective_value => JuMP.objective_value(env.master.model),
        :relative_gap => JuMP.relative_gap(env.master.model),
        :preprocessing_time => preprocessing_time,
        :solve_time => solve_time
    )

    if solution_strategy.verbose
        @info "Solver statistics" results...
    end

end