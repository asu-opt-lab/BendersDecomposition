function lazy_callback(cb_data, env::BendersEnv, cut_strategy::CutStrategy)
    status = JuMP.callback_node_status(cb_data, env.master.model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER    
        env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
        env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
        solve_sub!(env.sub, env.master.x_value)
        cuts, sub_obj_value = generate_cuts(env, cut_strategy)
        add_cuts!(env, cuts, sub_obj_value, cb_data)
    end
end

function user_callback(cb_data, env::BendersEnv, cut_strategy::CutStrategy)
    status = JuMP.callback_node_status(cb_data, env.master.model)
    depth = Ref{CPXLONG}()
    ret = CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEDEPTH, depth)
    n_count = Ref{CPXINT}()
    ret1 = CPXcallbackgetinfoint(cb_data, CPXCALLBACKINFO_NODECOUNT, n_count)

    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL 
        num_of_fraction_node += 1
        # println("cplex depth: $(depth[]), node count: $(n_count[]), num_of_fraction_node: $(num_of_fraction_node)")
        if num_of_fraction_node >= 1500
            num_of_fraction_node = 0
            explored_node = n_count[]
            env.master.x_value = JuMP.callback_value.(cb_data, env.master.var[:x])
            env.master.t_value = JuMP.callback_value.(cb_data, env.master.var[:t])
            println("cplex depth: $(depth[]), node count: $(n_count[])")

            zeros_indices = findall(x -> isapprox(x, 0.0; atol=1e-6), env.master.x_value)
            ones_indices = findall(x -> isapprox(x, 1.0; atol=1e-6), env.master.x_value)

            println("Indices where lb=ub=0: $(length(zeros_indices))")
            println("Indices where lb=ub=1: $(length(ones_indices))")

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