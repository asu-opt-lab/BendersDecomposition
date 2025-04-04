function generate_cuts(env::DisjunctiveBendersEnv, cut_strategy::DisjunctiveCut)

    sub_obj_val = get_subproblem_value(env.sub, env.master.integer_variable_values, cut_strategy.base_cut_strategy)
    # sub_obj_val = 1e06 # not use 

    disjunctive_inequality = select_disjunctive_inequality(env.master.integer_variable_values)

    update_dcglp!(env.dcglp, disjunctive_inequality, cut_strategy)

    solve_dcglp!(env, cut_strategy)
    
    cuts = merge_cuts(env, cut_strategy)

    println("--------------------------------DCGLP cut generation finished--------------------------------")

    return cuts, sub_obj_val
end

"""
    solve_dcglp!(env::AbstractBendersEnv, cut_strategy::DisjunctiveCut)

Solve the DCGLP problem and generate cuts.
"""
function solve_dcglp!(env::DisjunctiveBendersEnv, cut_strategy::DisjunctiveCut)
    # Initialize state and get values
    state = DCGLPState()
    x_value, t_value = env.master.integer_variable_values, env.master.continuous_variable_values
    
    # Update model with current values
    set_normalized_rhs.(env.dcglp.model[:conx], x_value)
    set_normalized_rhs.(env.dcglp.model[:cont], t_value)
    
    start_time = time()
    
    while true
        state.iteration += 1
        
        # Solve master problem and get values
        state.master_time += @elapsed begin
            k_values, v_values, other_values = solve_and_get_dcglp_values(env.dcglp.model, cut_strategy.norm_type)
        end
        
        # Solve subproblems and get cut coefficients
        state.sub_k_time += @elapsed begin
            dual_info_k, obj_value_k = solve_and_get_cut_coefficients(env.sub, k_values, cut_strategy.base_cut_strategy)
        end
        
        state.sub_v_time += @elapsed begin
            dual_info_v, obj_value_v = solve_and_get_cut_coefficients(env.sub, v_values, cut_strategy.base_cut_strategy)
        end
        
        state.total_time = time() - start_time
        
        # Update bounds and check termination
        update_bounds!(state, k_values, v_values, other_values, obj_value_k, obj_value_v, t_value, cut_strategy.norm_type)
        cut_strategy.verbose && print_dcglp_iteration_info(state)
        
        # Check termination conditions
        is_terminated(state, 30, 1000.0, 1e-6) && break
        
        # Add cuts if available
        !isempty(dual_info_k) && add_cuts!(env, dual_info_k, cut_strategy, K_CUTS)
        !isempty(dual_info_v) && add_cuts!(env, dual_info_v, cut_strategy, V_CUTS)
    end
end


function merge_cuts(env::DisjunctiveBendersEnv, cut_strategy::DisjunctiveCut)
    γ₀, γₓ, γₜ = generate_disjunctive_cuts(env.dcglp, cut_strategy.cut_strengthening_type)
    push!(env.dcglp.γ_values, (γ₀, γₓ, γₜ))
    master_disjunctive_cut = @expression(env.master.model, γ₀ - dot(γₓ, env.master.variables[:integer_variables]) - dot(γₜ, env.master.variables[:continuous_variables]))
    if cut_strategy.include_master_cuts
        push!(env.dcglp.master_cuts, master_disjunctive_cut)
        return env.dcglp.master_cuts
    end
    return master_disjunctive_cut
end
