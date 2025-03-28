function solve!(env::BendersEnv, ::Callback, cut_strategy::CutStrategy, params::BendersParams)

    start_time = time()
    df_root_node_preprocessing = root_node_preprocessing!(env, cut_strategy, params)
    params.time_limit -= df_root_node_preprocessing.total_time[end]
    number_of_subproblem_solves = 0

    function lazy_callback(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER
            number_of_subproblem_solves += 1   

            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])

            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_value = generate_cuts(env, cut_strategy)
            add_cuts!(env, cuts, sub_obj_value, cb_data)
        end
    end

    # Use the closure callbacks
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)
    set_optimizer_attribute(env.master.model, "CPX_PARAM_EPGAP", 1e-9)
    set_optimizer_attribute(env.master.model, "CPX_PARAM_EPAGAP", 0.0)
    set_optimizer_attribute(env.master.model, "CPX_PARAM_EPINT", 1e-7)

    MOI.set(env.master.model, MOI.RelativeGapTolerance(), 1e-9) 
    set_time_limit_sec(env.master.model, params.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)
    JuMP.optimize!(env.master.model)

    df_callback = DataFrame(
        node_count = JuMP.node_count(env.master.model),
        elapsed_time = time() - start_time,
        pure_callback_time = params.time_limit,
        objective_bound = JuMP.objective_bound(env.master.model),
        objective_value = JuMP.objective_value(env.master.model),
        relative_gap = JuMP.relative_gap(env.master.model),
        num_lazy = number_of_subproblem_solves,
        termination_status = termination_status(env.master.model)
    )
    
    return df_root_node_preprocessing, df_callback
end



function solve!(env::BendersEnv, ::Callback, cut_strategy::DisjunctiveCut, params::BendersParams)
    start_time = time()
    
    df_root_node_preprocessing = root_node_preprocessing!(env, cut_strategy.base_cut_strategy, params)
    params.time_limit -= df_root_node_preprocessing.total_time[end]

    # Dynamic information
    uc_spent_time, number_of_frac_solves = 0.0, 0

    # Static information
    number_of_subproblem_solves, number_of_user = 0, 0
    dcglp_times = []

    function lazy_callback(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER    
            number_of_subproblem_solves += 1

            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
            
            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_value = generate_cuts(env, cut_strategy.base_cut_strategy)
            add_cuts!(env, cuts, sub_obj_value, cb_data)
        end
    end

    function user_callback_lifting(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)

        if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            number_of_frac_solves += 1  

            if number_of_frac_solves >= 500
                # number_of_user += 1
                uc_start_time = time()

                depth = Ref{CPXLONG}()
                n_count = Ref{CPXINT}()
                @assert CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEDEPTH, depth) == 0
                @assert CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count) == 0

                println("depth: ", depth[], "node_count: ", n_count[])
                number_of_frac_solves = 0      

                env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
                env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])

                # Lifting 
                # Non-approximate
                # lb, ub = fill(NaN, length(env.master.var[:x])+length(env.master.var[:t])), fill(NaN, length(env.master.var[:x])+length(env.master.var[:t]))
                # @assert CPXcallbackgetlocallb(cb_data, lb, 0, length(lb) - 1) == 0
                # @assert CPXcallbackgetlocalub(cb_data, ub, 0, length(ub) - 1) == 0
                # zeros_indices = findall(i -> isapprox(lb[i], 0.0; atol=1e-6) && isapprox(ub[i], 0.0; atol=1e-6), 1:length(lb))
                # ones_indices = findall(i -> isapprox(lb[i], 1.0; atol=1e-6) && isapprox(ub[i], 1.0; atol=1e-6), 1:length(lb))
                
                # Approximate 
                zeros_indices = findall(x -> isapprox(x, 0.0; atol=1e-6), env.master.x_value)
                ones_indices = findall(x -> isapprox(x, 1.0; atol=1e-6), env.master.x_value)
                println("Indices where lb=ub=0: $(length(zeros_indices))")
                println("Indices where lb=ub=1: $(length(ones_indices))")
                
                solve_sub!(env.sub, env.master.x_value)
                dcglp_start_time = time()

                if length(zeros_indices) != 0 && length(ones_indices) != 0
                    cuts, sub_obj_value = generate_cuts_lifting(env, cut_strategy, zeros_indices, ones_indices)
                else
                    cuts, sub_obj_value = generate_cuts(env, cut_strategy)
                end

                # No lifitng
                # cuts, sub_obj_value = generate_cuts(env, cut_strategy)

                dcglp_end_time = time()
                dcglp_spent_time = dcglp_end_time - dcglp_start_time
                push!(dcglp_times, dcglp_spent_time)

                for _cut in cuts
                    cut = @build_constraint(0 >= _cut)
                    MOI.submit(env.master.model, MOI.UserCut(cb_data), cut)
                    number_of_user += 1
                end
                        
                uc_end_time = time()
                uc_spent_time += (uc_end_time - uc_start_time)
            end 
        end
    end

    function user_callback_lifting_local(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)

        if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            number_of_frac_solves += 1  

            if number_of_frac_solves >= 500
                number_of_user += 1
                uc_start_time = time()

                depth = Ref{CPXLONG}()
                n_count = Ref{CPXINT}()
                @assert CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEDEPTH, depth) == 0
                @assert CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count) == 0

                println("depth: ", depth[], "node_count: ", n_count[])
                number_of_frac_solves = 0      

                env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
                env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])

                # Lifting 
                # Non-approximate
                lb, ub = fill(NaN, length(env.master.var[:x])+length(env.master.var[:t])), fill(NaN, length(env.master.var[:x])+length(env.master.var[:t]))
                @assert CPXcallbackgetlocallb(cb_data, lb, 0, length(lb) - 1) == 0
                @assert CPXcallbackgetlocalub(cb_data, ub, 0, length(ub) - 1) == 0
                zeros_indices = findall(i -> isapprox(lb[i], 0.0; atol=1e-6) && isapprox(ub[i], 0.0; atol=1e-6), 1:length(lb))
                ones_indices = findall(i -> isapprox(lb[i], 1.0; atol=1e-6) && isapprox(ub[i], 1.0; atol=1e-6), 1:length(lb))
                
                # Approximate 
                # zeros_indices = findall(x -> isapprox(x, 0.0; atol=1e-6), env.master.x_value)
                # ones_indices = findall(x -> isapprox(x, 1.0; atol=1e-6), env.master.x_value)
                # println("Indices where lb=ub=0: $(length(zeros_indices))")
                # println("Indices where lb=ub=1: $(length(ones_indices))")
                
                solve_sub!(env.sub, env.master.x_value)
                dcglp_start_time = time()

                if length(zeros_indices) != 0 && length(ones_indices) != 0
                    cuts, sub_obj_value = generate_cuts_lifting(env, cut_strategy, zeros_indices, ones_indices)
                else
                    cuts, sub_obj_value = generate_cuts(env, cut_strategy)
                end

                # No lifitng
                # cuts, sub_obj_value = generate_cuts(env, cut_strategy)

                dcglp_end_time = time()
                dcglp_spent_time = dcglp_end_time - dcglp_start_time
                push!(dcglp_times, dcglp_spent_time)

                # disjunctive (global), all_found (local)
                n = 1
                if !isempty(cuts)
                    for _cut in cuts
                        cut = @build_constraint(0 >= _cut)
                        if n < length(cuts)
                            submit_local(env.master.model, MOI.UserCut(cb_data), cut)
                        else
                            @info _cut
                            MOI.submit(env.master.model, MOI.UserCut(cb_data), cut)
                        end
                        n += 1
                    end
                end

                uc_end_time = time()
                uc_spent_time += (uc_end_time - uc_start_time)
            end 
        end
    end
    
    # Use the closure callbacks
    set_binary.(env.master.model[:x])
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)
    set_attribute(env.master.model, MOI.UserCutCallback(), user_callback_lifting)
    # set_attribute(env.master.model, MOI.UserCutCallback(), user_callback_lifting_local)

    set_optimizer_attribute(env.master.model, "CPX_PARAM_EPGAP", 1e-9)
    set_optimizer_attribute(env.master.model, "CPX_PARAM_EPAGAP", 0.0)
    set_optimizer_attribute(env.master.model, "CPX_PARAM_EPINT", 1e-7)

    MOI.set(env.master.model, MOI.RelativeGapTolerance(), 1e-9) 
    set_time_limit_sec(env.master.model, params.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)
    set_optimizer_attribute(env.master.model, "CPXPARAM_MIP_Display", 3)
    set_optimizer_attribute(env.master.model, "CPX_PARAM_RANDOMSEED", 1218)

    JuMP.optimize!(env.master.model)
    println("number of user cuts $number_of_user")
    df_callback = DataFrame(
        node_count = JuMP.node_count(env.master.model),
        elapsed_time = time() - start_time,
        pure_callback_time = params.time_limit,
        objective_bound = JuMP.objective_bound(env.master.model),
        objective_value = JuMP.objective_value(env.master.model),
        relative_gap = JuMP.relative_gap(env.master.model),
        num_lazy = number_of_subproblem_solves,
        average_dcglp_time = isempty(dcglp_times) ? -1 : sum(dcglp_times)/length(dcglp_times),
        strengthen_percentage = isempty(env.dcglp.strengthen_used) ? -1 : sum(env.dcglp.strengthen_used)/length(env.dcglp.strengthen_used),
        user_cut_spent_time = uc_spent_time,
        num_user_cut = number_of_user,
        termination_status = termination_status(env.master.model)
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