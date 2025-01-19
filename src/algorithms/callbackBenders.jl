function solve!(env::BendersEnv, ::Callback, cut_strategy::CutStrategy, params::BendersParams)

    start_time = time()
    # time_limit = params.time_limit
    # params.time_limit = 600
    # df_root_node_preprocessing = root_node_preprocessing!(env, cut_strategy, params)
    # params.time_limit = time_limit
    # params.time_limit -= df_root_node_preprocessing.total_time[end]

    function lazy_callback(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER         
            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
            
            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_value = generate_cuts(env, cut_strategy)
            add_cuts!(env, cuts, sub_obj_value, cb_data)
        end
    end

    # Use the closure callbacks
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)
    MOI.set(env.master.model, MOI.RelativeGapTolerance(), params.gap_tolerance*1e-2) # make sure the relative gap is small
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



function solve!(env::BendersEnv, ::Callback, cut_strategy::DisjunctiveCut, params::BendersParams)
    start_time = time()
    # time_limit = params.time_limit
    # params.time_limit = 600
    # df_root_node_preprocessing = root_node_preprocessing!(env, cut_strategy.base_cut_strategy, params)
    # params.time_limit = time_limit
    # params.time_limit -= df_root_node_preprocessing.total_time[end]
    # df2 = root_node_preprocessing!(env, cut_strategy, params)
    number_of_subproblem_solves = 0

    function lazy_callback(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER     

            number_of_subproblem_solves += 1

            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            value_t = JuMP.callback_value.(cb_data, env.master.var[:t])
            
            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_value = generate_cuts(env, cut_strategy.base_cut_strategy)
            if value_t <= sub_obj_value - 1e-06
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
        if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL && depth[] <= 10 && number_of_subproblem_solves >= 200 
            number_of_subproblem_solves = 0      
            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
            
            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_value = generate_cuts(env, cut_strategy)
            if env.master.t_value <= sub_obj_value - 1e-06
                for _cut in cuts
                    cut = @build_constraint(0 .>= _cut)
                    # @info "User cut" cut
                    MOI.submit.(env.master.model, MOI.UserCut(cb_data), cut)
                end
            end
        end
    end

    # Use the closure callbacks
    set_binary.(env.master.model[:x])
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)
    set_attribute(env.master.model, MOI.UserCutCallback(), user_callback)

    MOI.set(env.master.model, MOI.RelativeGapTolerance(), params.gap_tolerance*1e-2) # make sure the relative gap is small
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




# for base cut strategy
function root_node_preprocessing!(env::BendersEnv, cut_strategy::CutStrategy, params::BendersParams)
    relax_integrality(env.master.model)
    df = solve!(env, Sequential(), cut_strategy, params)
    set_binary.(env.master.model[:x])
    return df
end

function add_cuts!(env::BendersEnv, expressions::Vector{Any}, sub_obj_values::Vector{Float64}, cb_data)
    for (idx, (expr, sub_obj)) in enumerate(zip(expressions, sub_obj_values))
        if env.master.t_value[idx] <= sub_obj - 1e-06
            cut = @build_constraint(0 >= expr)
            MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut)
        end
    end
end

function add_cuts!(env::BendersEnv, expression::Any, sub_obj_value::Float64, cb_data)
    if env.master.t_value <= sub_obj_value - 1e-06
        cut = @build_constraint(0 >= expression)
        MOI.submit(env.master.model, MOI.LazyConstraint(cb_data), cut)
    end
end
