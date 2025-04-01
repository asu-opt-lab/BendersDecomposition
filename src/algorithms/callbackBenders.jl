function solve!(env::BendersEnv, ::Callback, cut_strategy::CutStrategy, params::BendersParams)

    start_time = time()
    time_limit = params.time_limit
    params.time_limit = 600
    df_root_node_preprocessing = root_node_preprocessing!(env, cut_strategy, params)
    params.time_limit = time_limit
    params.time_limit -= df_root_node_preprocessing.total_time[end]

    function lazy_callback(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER 
            n_count = Ref{CPXINT}()
            ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)

            # @info "lazy constraints node: $(n_count[])"        
            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
            
            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_value = generate_cuts(env, cut_strategy)
            add_cuts!(env, cuts, sub_obj_value, cb_data)
        end
    end

    # Use the closure callbacks
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)
    # MOI.set(env.master.model, MOI.RelativeGapTolerance(), params.gap_tolerance*1e-2) # make sure the relative gap is small
    set_time_limit_sec(env.master.model, params.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)
    set_optimizer_attribute(env.master.model, "CPXPARAM_MIP_Display", 3)

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
    time_limit = params.time_limit
    params.time_limit = 600
    df_root_node_preprocessing = root_node_preprocessing!(env, cut_strategy.base_cut_strategy, params)
    params.time_limit = time_limit
    params.time_limit -= df_root_node_preprocessing.total_time[end]
    # df2 = root_node_preprocessing!(env, cut_strategy, params)
    number_of_subproblem_solves = 0
    explored_node = 0
    cut_count = 0
    num_of_fraction_node = 0

    function lazy_callback(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        if status == MOI.CALLBACK_NODE_STATUS_INTEGER     
            n_count = Ref{CPXINT}()
            ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)
            number_of_subproblem_solves += 1
            # @info "lazy constraints node: $(n_count[])"
            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
            set_optimizer_attribute(env.sub.model, "CPX_PARAM_LPMETHOD", 0)
            set_optimizer_attribute(env.sub.model, "CPX_PARAM_EPOPT", 1e-06)

            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_value = generate_cuts(env, cut_strategy.base_cut_strategy)
            add_cuts!(env, cuts, sub_obj_value, cb_data)
        end
    end

    # function user_callback(cb_data)
    #     status = JuMP.callback_node_status(cb_data, env.master.model)
    #     depth = Ref{CPXLONG}()
    #     if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL && depth[] <= 10 && number_of_subproblem_solves >= 200 
    #         number_of_subproblem_solves = 0      
    #         env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
    #         env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
            
    #         # solve_sub!(env.sub, env.master.x_value)
    #         cuts, sub_obj_value = generate_cuts(env, cut_strategy)
    #         for _cut in cuts
    #             cut = @build_constraint(0 >= _cut)
    #             MOI.submit(env.master.model, MOI.UserCut(cb_data), cut)
    #         end

    #     end
    # end

    # function user_callback_lifting(cb_data)
    #     status = JuMP.callback_node_status(cb_data, env.master.model)
    #     depth = Ref{CPXLONG}()
    #     ret = CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEDEPTH, depth)
    #     n_count = Ref{CPXINT}()
    #     ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)

    #     if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL && 
    #         10 <= depth[] <= 30 && 
    #         number_of_subproblem_solves >= 20

    #         number_of_subproblem_solves = 0
    #         explored_node = n_count[]
    #         env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
    #         env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
    #         lb, ub = fill(NaN, length(env.master.var[:x])+length(env.master.var[:t])), fill(NaN, length(env.master.var[:x])+length(env.master.var[:t]))
    #         @assert CPXcallbackgetlocallb(cb_data, lb, 0, length(lb) - 1) == 0
    #         @assert CPXcallbackgetlocalub(cb_data, ub, 0, length(ub) - 1) == 0
    #         # @info "lb" lb
    #         # @info "ub" ub
    #         println("cplex depth: $(depth[])")
    #         zeros_indices = findall(i -> isapprox(lb[i], 0.0; atol=1e-6) && isapprox(ub[i], 0.0; atol=1e-6), 1:length(lb))           
    #         ones_indices = findall(i -> isapprox(lb[i], 1.0; atol=1e-6) && isapprox(ub[i], 1.0; atol=1e-6), 1:length(lb))
            
    #         println("Indices where lb=ub=0: $(length(zeros_indices))")
    #         println("Indices where lb=ub=1: $(length(ones_indices))")


    #         solve_sub!(env.sub, env.master.x_value)
    #         cuts, sub_obj_value = generate_cuts_lifting(env, cut_strategy, zeros_indices, ones_indices)

    #         for _cut in cuts
    #             cut = @build_constraint(0 >= _cut)
    #             MOI.submit(env.master.model, MOI.UserCut(cb_data), cut)
    #         end
    #     end
    # end

    # function user_callback_lifting_1(cb_data)
    #     status = JuMP.callback_node_status(cb_data, env.master.model)
    #     depth = Ref{CPXLONG}()
    #     ret = CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEDEPTH, depth)
    #     n_count = Ref{CPXINT}()
    #     ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)

    #     if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL 
    #         num_of_fraction_node += 1
    #         if num_of_fraction_node >= 2500
    #             num_of_fraction_node = 0
    #             explored_node = n_count[]
    #             env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
    #             env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
    #             lb, ub = fill(NaN, length(env.master.var[:x])+length(env.master.var[:t])), fill(NaN, length(env.master.var[:x])+length(env.master.var[:t]))
    #             @assert CPXcallbackgetlocallb(cb_data, lb, 0, length(lb) - 1) == 0
    #             @assert CPXcallbackgetlocalub(cb_data, ub, 0, length(ub) - 1) == 0
    #             # @info "lb" lb
    #             # @info "ub" ub
    #             println("cplex depth: $(depth[])")
    #             zeros_indices = findall(i -> isapprox(lb[i], 0.0; atol=1e-6) && isapprox(ub[i], 0.0; atol=1e-6), 1:length(lb))           
    #             ones_indices = findall(i -> isapprox(lb[i], 1.0; atol=1e-6) && isapprox(ub[i], 1.0; atol=1e-6), 1:length(lb))
                
    #             println("Indices where lb=ub=0: $(length(zeros_indices))")
    #             println("Indices where lb=ub=1: $(length(ones_indices))")


    #             solve_sub!(env.sub, env.master.x_value)
    #             if zeros_indices == [] || ones_indices == []
    #                 cuts, sub_obj_value = generate_cuts(env, cut_strategy)
    #             else
    #                 cuts, sub_obj_value = generate_cuts_lifting(env, cut_strategy, zeros_indices, ones_indices)
    #                 # cuts, sub_obj_value = generate_cuts(env, cut_strategy)
    #             end

    #             for _cut in cuts
    #                 cut = @build_constraint(0 >= _cut)
    #                 MOI.submit(env.master.model, MOI.UserCut(cb_data), cut)
    #             end
    #         end
    #     end
    # end

    function user_callback_lifting_approx(cb_data)
        status = JuMP.callback_node_status(cb_data, env.master.model)
        depth = Ref{CPXLONG}()
        ret = CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEDEPTH, depth)
        n_count = Ref{CPXINT}()
        ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)

        if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL 
            num_of_fraction_node += 1
            # println("cplex depth: $(depth[]), node count: $(n_count[]), num_of_fraction_node: $(num_of_fraction_node)")
            if num_of_fraction_node >= 100
                num_of_fraction_node = 0
                explored_node = n_count[]
                env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
                env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
                println("cplex depth: $(depth[]), node count: $(n_count[])")

                zeros_indices = findall(x -> isapprox(x, 0.0; atol=1e-6), env.master.x_value)
                ones_indices = findall(x -> isapprox(x, 1.0; atol=1e-6), env.master.x_value)

                println("Indices where lb=ub=0: $(length(zeros_indices))")
                println("Indices where lb=ub=1: $(length(ones_indices))")
                
                set_optimizer_attribute(env.sub.model, "CPX_PARAM_LPMETHOD", 0)
                set_optimizer_attribute(env.sub.model, "CPX_PARAM_EPOPT", 1e-06)
                set_optimizer_attribute(env.sub.model, "CPX_PARAM_ITLIM", 9223372036800000000)
                set_optimizer_attribute(env.sub.model, MOI.Silent(), true) 
                # set_optimizer_attribute(env.sub.model, "CPXPARAM_MIP_Display", 3)    
                solve_sub!(env.sub, env.master.x_value)
                if zeros_indices == [] || ones_indices == []
                    cuts, sub_obj_value = generate_cuts(env, cut_strategy)
                else
                    cuts, sub_obj_value = generate_cuts_lifting(env, cut_strategy, zeros_indices, ones_indices)
                end
                # cuts, sub_obj_value = generate_cuts(env, cut_strategy)
                # cuts, sub_obj_value = generate_cuts_lifting(env, cut_strategy, zeros_indices, ones_indices)
                for _cut in cuts
                    cut = @build_constraint(0 >= _cut)
                    MOI.submit(env.master.model, MOI.UserCut(cb_data), cut)
                end
            end
        end
    end
    # Use the closure callbacks
    set_binary.(env.master.model[:x])
    set_attribute(env.master.model, MOI.LazyConstraintCallback(), lazy_callback)
    # set_attribute(env.master.model, MOI.UserCutCallback(), user_callback)
    # set_attribute(env.master.model, MOI.UserCutCallback(), user_callback_lifting_1)
    set_attribute(env.master.model, MOI.UserCutCallback(), user_callback_lifting_approx)

    # MOI.set(env.master.model, MOI.RelativeGapTolerance(), params.gap_tolerance*1e-2) # make sure the relative gap is small
    set_time_limit_sec(env.master.model, params.time_limit)
    set_optimizer_attribute(env.master.model, MOI.Silent(), false)
    set_optimizer_attribute(env.master.model, "CPXPARAM_MIP_Display", 3)
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_EPINT", 0.0)
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_EPGAP", 1e-9)
    # set_optimizer_attribute(env.master.model, "CPX_PARAM_EPRHS", 1e-9)

    # @info "number of all constraints" length(all_constraints(env.master.model, include_variable_in_set_constraints=false))
    JuMP.optimize!(env.master.model)
    @info "node count" JuMP.node_count(env.master.model)
    @info "elapsed time" time() - start_time
    @info "objective bound" JuMP.objective_bound(env.master.model)
    @info "objective value" JuMP.objective_value(env.master.model)
    @info "relative gap" JuMP.relative_gap(env.master.model)
    # @info "number of all constraints" length(all_constraints(env.master.model, include_variable_in_set_constraints=false))
    return JuMP.objective_value(env.master.model),time() - start_time
end




# for base cut strategy
function root_node_preprocessing!(env::BendersEnv, cut_strategy::CutStrategy, params::BendersParams)
    relax_integrality(env.master.model)
    df = solve!(env, Sequential(), cut_strategy, params)
    set_binary.(env.master.model[:x])
    println("--------------------------------Root node preprocessing finished--------------------------------")
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
