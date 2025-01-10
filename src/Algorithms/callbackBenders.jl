function solve!(env::BendersEnv, ::Callback, cut_strategy::CutStrategy, params::BendersParams)

    start_time = time()
    df_root_node_preprocessing = root_node_preprocessing!(env, cut_strategy, params)
    # params.time_limit = 7200.0 # set time for callback (internal parameter)
    params.time_limit -= df_root_node_preprocessing.total_time[end]

    function lazy_callback(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER         
            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            value_t = JuMP.callback_value.(cb_data, env.master.var[:t])
            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_value = generate_cuts(env, cut_strategy)
            if sum(value_t) <= sub_obj_value - params.gap_tolerance
                for _cut in cuts
                    cut = @build_constraint(0 >= _cut)
                    MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut)
                end
            end
        end
    end

    # Use the closure callbacks
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)
    # MOI.set(env.master.model, MOI.RelativeGapTolerance(), params.gap_tolerance) # Not necessary if applied in config file
    set_time_limit_sec(env.master.model, params.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)

    JuMP.optimize!(env.master.model)

    @info "node count" JuMP.node_count(env.master.model)
    @info "elapsed time" time() - start_time
    @info "pure callback time", params.time_limit
    @info "objective bound" JuMP.objective_bound(env.master.model)
    @info "objective value" JuMP.objective_value(env.master.model)
    @info "relative gap" JuMP.relative_gap(env.master.model)

    df_callback = DataFrame(
        node_count = JuMP.node_count(env.master.model),
        elapsed_time = time() - start_time,
        pure_callback_time = params.time_limit,
        objective_bound = JuMP.objective_bound(env.master.model),
        objective_value = JuMP.objective_value(env.master.model),
        relative_gap = JuMP.relative_gap(env.master.model)
    )
    
    return df_root_node_preprocessing, df_callback
end

function solve!(env::BendersEnv, ::Callback, cut_strategy::DisjunctiveCut, params::BendersParams)
    start_time = time()
    df_root_node_preprocessing = root_node_preprocessing!(env, cut_strategy, params)
    # params.time_limit = 7200.0 # set time for callback (internal parameter)
    params.time_limit -= df_root_node_preprocessing.total_time[end]
    number_of_subproblem_solves = 0

    function lazy_callback(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER     

            number_of_subproblem_solves += 1

            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            value_t = JuMP.callback_value.(cb_data, env.master.var[:t])
            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_value = generate_cuts(env, cut_strategy.base_cut_strategy)
            
            if sum(value_t) <= sub_obj_value - params.gap_tolerance
                for _cut in cuts
                    cut = @build_constraint(0 >= _cut)
                    MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut)
                end
            end
        end
    end

    function user_callback(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        depth = Ref{CPXLONG}()
        ret = CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEDEPTH, depth)
        
        n_count = Ref{CPXINT}()
        ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)

        if ret == 0 && ret1 == 0
            if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL 
                if depth[] <= 5 
                # if depth[] <= 5 || (depth[] > 5 && n_count[] % 100 == 0 && n_count[] != 0)
                    @info "depth", depth[]
                    @info "node_count", n_count[]
                    number_of_subproblem_solves = 0      
                    env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
                    env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
                    solve_sub!(env.sub, env.master.x_value)
                    cuts, sub_obj_value = generate_cuts(env, cut_strategy)
                    if sum(env.master.t_value) <= sub_obj_value - params.gap_tolerance
                        for _cut in cuts
                            cut = @build_constraint(0 .>= _cut)
                            # @info "User cut" cut
                            MOI.submit.(env.master.model, MOI.UserCut(cb_data), cut)
                        end
                    end
                end
            end
        else
            @info "Wrong value for ret or ret1"
            exit(1)
        end
    end

    # Use the closure callbacks
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)
    set_attribute(env.master.model, MOI.UserCutCallback(), user_callback)

    # MOI.set(env.master.model, MOI.RelativeGapTolerance(), params.gap_tolerance) # Not necessary if applied in config file
    set_time_limit_sec(env.master.model, params.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)


    JuMP.optimize!(env.master.model)
    @info "node count" JuMP.node_count(env.master.model)
    @info "elapsed time" time() - start_time
    @info "pure callback time", params.time_limit
    @info "objective bound" JuMP.objective_bound(env.master.model)
    @info "objective value" JuMP.objective_value(env.master.model)
    @info "relative gap" JuMP.relative_gap(env.master.model)

    df_callback = DataFrame(
        node_count = JuMP.node_count(env.master.model),
        elapsed_time = time() - start_time,
        pure_callback_time = params.time_limit,
        objective_bound = JuMP.objective_bound(env.master.model),
        objective_value = JuMP.objective_value(env.master.model),
        relative_gap = JuMP.relative_gap(env.master.model)
    )
    
    return df_root_node_preprocessing, df_callback
end


# for base cut strategy
function root_node_preprocessing!(env::BendersEnv, cut_strategy::CutStrategy, params::BendersParams)
    relax_integrality(env.master.model)
    df = solve!(env, Sequential(), cut_strategy, params)
    set_binary.(env.master.model[:x])
    return df
end