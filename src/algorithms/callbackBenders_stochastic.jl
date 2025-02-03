function solve!(env::BendersEnv, ::StochasticCallback, cut_strategy::CutStrategy, params::BendersParams)

    start_time = time()
    # time_limit = params.time_limit
    # params.time_limit = 600
    # df_root_node_preprocessing = root_node_stochastic_preprocessing!(env, cut_strategy, params)
    # params.time_limit = time_limit
    # params.time_limit -= df_root_node_preprocessing.total_time[end]

    function lazy_callback(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER         
            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            value_t = JuMP.callback_value.(cb_data, env.master.var[:t])
            
            for scenario in 1:env.data.num_scenarios
                solve_sub!(env.sub.sub_problems[scenario], env.master.x_value)
                cuts, sub_obj_value = generate_cuts(env, cut_strategy, scenario)
                if value_t[scenario] <= sub_obj_value - 1e-06
                    for _cut in cuts
                    cut = @build_constraint(0 >= _cut)
                        MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut)
                    end
                end
            end
        end
    end

    # Use the closure callbacks
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)
    MOI.set(env.master.model, MOI.RelativeGapTolerance(), params.gap_tolerance*1e-2) # make sure the relative gap is small

    # turn off CPLEX cuts
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_CUTSFACTOR", 0)
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_CLIQUES", -1)
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_COVERS", -1)
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_FLOWCOVERS", -1)
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_FRACCUTS", -1)
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_GUBCOVERS", -1)
    # # set_optimizer_attribute(env.master.model, "CPX_PARAM_IMPLIED", -1)
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_MCFCUTS", -1)
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_MIRCUTS", -1)
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_ZEROHALFCUTS", -1)

    set_time_limit_sec(env.master.model, params.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)
    JuMP.optimize!(env.master.model)
    @info "node count" JuMP.node_count(env.master.model)
    @info "elapsed time" time() - start_time
    @info "objective bound" JuMP.objective_bound(env.master.model)
    @info "objective value" JuMP.objective_value(env.master.model)
    @info "relative gap" JuMP.relative_gap(env.master.model)
    return JuMP.objective_value(env.master.model),time() - start_time
end



function root_node_stochastic_preprocessing!(env::BendersEnv, cut_strategy::CutStrategy, params::BendersParams)
    relax_integrality(env.master.model)
    df = solve!(env, StochasticSequential(), cut_strategy, params)
    set_binary.(env.master.model[:x])
    return df
end