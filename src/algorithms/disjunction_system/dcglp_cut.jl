function generate_cuts(env::BendersEnv, cut_strategy::DisjunctiveCut)

    sub_obj_val = get_subproblem_value(env) 

    disjunctive_inequality = select_disjunctive_inequality(env.master.x_value)

    update_dcglp!(env.dcglp, disjunctive_inequality, cut_strategy)
    
    solve_dcglp!(env, cut_strategy)
    
    cuts = merge_cuts(env, cut_strategy)

    return cuts, sub_obj_val
end

function solve_dcglp!(env::BendersEnv, cut_strategy::DisjunctiveCut)

    log = DCGLPIterationLog()
    state = DCGLPState()

    x_value, t_value = env.master.x_value, env.master.t_value

    set_normalized_rhs.(env.dcglp.model[:conx], x_value)
    set_normalized_rhs.(env.dcglp.model[:cont], t_value)
    log.start_time = time()

    while true

        state.iteration += 1

        master_time = @elapsed begin
            k_values, v_values, other_values = solve_and_get_dcglp_values(env.dcglp.model, cut_strategy.norm_type)
        end
        log.master_time += master_time

        sub_k_time = @elapsed begin
            if_add_cuts_k, dual_info_k, obj_value_k = generate_cut_coefficients(env.sub, k_values, cut_strategy.base_cut_strategy)
        end
        log.sub_k_time += sub_k_time
        
        sub_v_time = @elapsed begin
            if_add_cuts_v, dual_info_v, obj_value_v = generate_cut_coefficients(env.sub, v_values, cut_strategy.base_cut_strategy)
        end
        log.sub_v_time += sub_v_time

        update_bounds!(state, k_values, v_values, other_values, obj_value_k, obj_value_v, t_value, cut_strategy.norm_type)
        
        record_iteration!(log, state)

        cut_strategy.verbose && print_dcglp_iteration_info(state, log)

        is_terminated(state, log) && break
        
        if if_add_cuts_k
            add_cuts_k!(env, dual_info_k, cut_strategy)
        end
        if if_add_cuts_v
            add_cuts_v!(env, dual_info_v, cut_strategy)
        end
        
    end

end


function is_terminated(state, log)
    return state.gap <= 1e-3  || state.UB - state.LB <= 1e-03 || (state.UB_k <= 1e-3 && state.UB_v <= 1e-3) || get_total_time(log) >= 200 || state.iteration >= 50
end


function print_dcglp_iteration_info(state, log)
    @printf("   Iter: %4d | LB: %12.4f | UB: %11.4f | Gap: %8.2f%% | UB_k: %11.4f | UB_v: %11.4f \n",
           state.iteration, state.LB, state.UB, state.gap, state.UB_k, state.UB_v)
end

function update_bounds!(state, k_values, v_values, other_values, obj_value_k, obj_value_v, t_value, norm_type::LNorm)
    state.LB = other_values.τ
    diff_st = obj_value_k .+ obj_value_v .- t_value
    state.UB = update_UB!(state.UB, other_values.sx, diff_st, norm_type)
    state.UB_k = obj_value_k .- k_values.t
    state.UB_v = obj_value_v .- v_values.t
    state.gap = (state.UB - state.LB)/abs(state.UB) * 100
end

update_UB!(UB,_sx,diff_st,::L1Norm) = min(UB,norm([ _sx; diff_st], Inf))
update_UB!(UB,_sx,diff_st,::L2Norm) = min(UB,norm([ _sx; diff_st], 2))
update_UB!(UB,_sx,diff_st,::LInfNorm) = min(UB,norm([ _sx; diff_st], 1))

function merge_cuts(env::BendersEnv, cut_strategy::DisjunctiveCut)
    γ₀, γₓ, γₜ = generate_disjunctive_cuts(env.dcglp, cut_strategy.cut_strengthening_type)
    push!(env.dcglp.γ_values, (γ₀, γₓ, γₜ))
    master_disjunctive_cut = @expression(env.master.model, γ₀ + dot(γₓ, env.master.var[:x]) + dot(γₜ, env.master.var[:t]))
    if cut_strategy.include_master_cuts
        push!(env.dcglp.master_cuts, [master_disjunctive_cut])
        return env.dcglp.master_cuts
    end
    return master_disjunctive_cut
end

# ============================================================================
# Generate strengthened cuts for the DCGLP
# ============================================================================

function generate_disjunctive_cuts(dcglp::DCGLP, ::PureDisjunctiveCut)
    optimize!(dcglp.model)
    γₜ = dual.(dcglp.model[:cont])
    γ₀ = dual(dcglp.model[:con0])
    γₓ = dual.(dcglp.model[:conx])
    return γ₀, γₓ, γₜ
end


function generate_disjunctive_cuts(dcglp::DCGLP, ::StrengthenedDisjunctiveCut)
    optimize!(dcglp.model)
    
    σ₁::Float64 = dual(dcglp.disjunctive_inequalities_constraints[1])
    σ₂::Float64 = dual(dcglp.disjunctive_inequalities_constraints[2])
    γₜ::Float64 = dual(dcglp.γ_constraints[:γₜ])
    γ₀::Float64 = dual(dcglp.γ_constraints[:γ₀])
    γₓ::Vector{Float64} = dual.(dcglp.γ_constraints[:γₓ])
    if iszero(σ₁) && iszero(σ₂)
        return γ₀, γₓ, γₜ
    end
    γ₁ = γₓ .- dual.(dcglp.model[:conv1])
    γ₂ = γₓ .- dual.(dcglp.model[:conv2])
    σ_sum = σ₂ + σ₁
    if !iszero(σ_sum)
        m = (γ₁ .- γ₂) / σ_sum
        m_lb = floor.(m)
        m_ub = ceil.(m)
        γₓ = min.(γ₁-σ₁*m_lb, γ₂+σ₂*m_ub)
    end
    
    return γ₀, γₓ, γₜ
end
