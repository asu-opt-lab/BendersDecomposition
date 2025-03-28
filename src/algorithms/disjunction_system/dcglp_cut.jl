function generate_cuts(env::BendersEnv, cut_strategy::DisjunctiveCut)

    sub_obj_val = get_subproblem_value(env.sub, env.master.x_value, cut_strategy.base_cut_strategy) 

    disjunctive_inequality = select_disjunctive_inequality(env.master.x_value)

    split_index = findall(x -> x == 1 , disjunctive_inequality[1])[1]
    push!(env.dcglp.split_index, split_index)

    update_dcglp!(env.dcglp, disjunctive_inequality, cut_strategy)
    
    solve_dcglp!(env, cut_strategy)
    
    cuts = merge_cuts(env, cut_strategy)
    println("--------------------------------DCGLP cut generation finished--------------------------------")
    return cuts, sub_obj_val
end

function solve_dcglp!(env::BendersEnv, cut_strategy::DisjunctiveCut)

    log = DCGLPIterationLog()
    state = DCGLPState()

    x_value, t_value = env.master.x_value, env.master.t_value

    set_normalized_rhs.(env.dcglp.model[:conx], x_value)
    set_normalized_rhs.(env.dcglp.model[:cont], t_value)
    log.start_time = time()

    # set time limitation 
    set_time_limit_sec(env.dcglp.model, 30.0)
    add_norm_cons(env)

    while true

        state.iteration += 1

        # println("before time: $(time() - log.start_time)")
        master_time = @elapsed begin
            k_values, v_values, other_values = solve_and_get_dcglp_values(env.dcglp.model, cut_strategy.norm_type)
        end
        log.master_time += master_time
        # println("after time: $(time() - log.start_time)")

        # set time limitation 
        if time() - log.start_time >= 30 || k_values == nothing
            break
        end

        sub_k_time = @elapsed begin
            dual_info_k, obj_value_k = generate_cut_coefficients(env.sub, k_values, cut_strategy.base_cut_strategy)
        end
        log.sub_k_time += sub_k_time
        # println("dual_info_k: $dual_info_k")
        sub_v_time = @elapsed begin
            dual_info_v, obj_value_v = generate_cut_coefficients(env.sub, v_values, cut_strategy.base_cut_strategy)
        end
        log.sub_v_time += sub_v_time
        # println("dual_info_v: $dual_info_v")
        update_bounds!(state, k_values, v_values, other_values, obj_value_k, obj_value_v, t_value, cut_strategy.norm_type)
        
        record_iteration!(log, state)

        cut_strategy.verbose && print_dcglp_iteration_info(state, log)

        is_terminated(state, log) && break
        
        if dual_info_k != []
            add_cuts_k!(env, dual_info_k, cut_strategy)
        end
        if dual_info_v != []
            add_cuts_v!(env, dual_info_v, cut_strategy)
        end
    end
    # k_x, v_x = value.(env.dcglp.var_kv_x[:kₓ]), value.(env.dcglp.var_kv_x[:vₓ])
    # println("k_x: $k_x", "v_x: $v_x")
    # println("k_x_index:", k_x[env.dcglp.split_index[end]], "v_x_index: ", v_x[env.dcglp.split_index[end]])
    # println("sum kappa_x and nu_x for disjunctive index", k_x[env.dcglp.split_index[end]] + v_x[env.dcglp.split_index[end]])
    # statistics_of_disjunctive_cuts(env)
end

# Stopping criteria: consecutive improvement
# function is_terminated(state, log)
#     if state.iteration == 1 
#         if state.gap <= 1e-3  || state.UB - state.LB <= 1e-03 || get_total_time(log) >= 30
#             return true
#         else
#             return false
#         end
#     else
#         if state.gap <= 1e-3  || state.UB - state.LB <= 1e-03 || get_total_time(log) >= 30 || (state.LB_set[end] - state.LB_set[end-1]) / state.LB_set[end-1] <= 0.1
#             return true
#         else        
#             return false
#         end
#     end
# end

# Stopping criteria: original
function is_terminated(state, log)
    return state.gap <= 1e-3  || state.UB - state.LB <= 1e-03 || get_total_time(log) >= 30 || state.iteration >= 30 #|| (state.UB_k <= 1e-3 && state.UB_v <= 1e-3) 
end

function print_dcglp_iteration_info(state, log)
    @printf("   Iter: %4d | LB: %8.4f | UB: %8.4f | Gap: %6.2f%% | UB_k: %8.2f | UB_v: %8.2f | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
           state.iteration, state.LB, state.UB, state.gap, sum(state.UB_k), sum(state.UB_v), log.master_time, log.sub_k_time, log.sub_v_time)
end

function update_bounds!(state, k_values, v_values, other_values, obj_value_k, obj_value_v, t_value, norm_type::LNorm)
    state.LB = other_values.τ 
    diff_st = (obj_value_k .+ obj_value_v .- t_value)
    state.UB = update_UB!(state.UB, other_values.sx, diff_st, norm_type)
    state.UB_k = obj_value_k .- k_values.t
    state.UB_v = obj_value_v .- v_values.t
    state.gap = (state.UB - state.LB)/abs(state.UB) * 100
    # @info "k_values", k_values
    # @info "v_values", v_values
    # @info "obj_value_k", obj_value_k
    # @info "obj_value_v", obj_value_v
    # @info "diff_st", diff_st
    push!(state.LB_set, other_values.τ)
end

update_UB!(UB,_sx,diff_st,::L1Norm) = min(UB,norm([ _sx; diff_st], Inf))
update_UB!(UB,_sx,diff_st,::L2Norm) = min(UB,norm([ _sx; diff_st], 2))
update_UB!(UB,_sx,diff_st,::LInfNorm) = min(UB,norm([ _sx; diff_st], 1))

function merge_cuts(env::BendersEnv, cut_strategy::DisjunctiveCut)
    if termination_status(env.dcglp.model) != OPTIMAL
        println("no disjunctive cuts generated")
        return env.dcglp.master_cuts
    end
    γ₀, γₓ, γₜ = generate_disjunctive_cuts(env.dcglp, cut_strategy.cut_strengthening_type)
    push!(env.dcglp.γ_values, (γ₀, γₓ, γₜ))
    master_disjunctive_cut = @expression(env.master.model, γ₀ + dot(γₓ, env.master.var[:x]) + dot(γₜ, env.master.var[:t]))
    if cut_strategy.include_master_cuts
        push!(env.dcglp.master_cuts, master_disjunctive_cut)
        return env.dcglp.master_cuts
    end
    return master_disjunctive_cut
end

# ============================================================================
# Generate strengthened cuts for the DCGLP
# ============================================================================

function generate_disjunctive_cuts(dcglp::DCGLP, ::PureDisjunctiveCut)
    # optimize!(dcglp.model)
    γₜ = dual.(dcglp.model[:cont])
    γ₀ = dual(dcglp.model[:con0])
    γₓ = dual.(dcglp.model[:conx])
    return γ₀, γₓ, γₜ
end


function generate_disjunctive_cuts(dcglp::DCGLP, ::StrengthenedDisjunctiveCut)
    # optimize!(dcglp.model)
    
    σ₁::Float64 = dual(dcglp.disjunctive_inequalities_constraints[1])
    σ₂::Float64 = dual(dcglp.disjunctive_inequalities_constraints[2])
    γₜ::Vector{Float64} = dual.(dcglp.γ_constraints[:γₜ])
    γ₀::Float64 = dual(dcglp.γ_constraints[:γ₀])
    γₓ::Vector{Float64} = -dual.(dcglp.γ_constraints[:γₓ])
    println("DCGLP Sigma Values: [σ₁: $σ₁, σ₂: $σ₂]")
    if abs(σ₁) <= 1e-6 && abs(σ₂) <= 1e-6
        push!(dcglp.strengthen_used, 0)
        return γ₀, -γₓ, γₜ
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
    push!(dcglp.strengthen_used, 1)
    return γ₀, -γₓ, γₜ
end

function statistics_of_disjunctive_cuts(env::BendersEnv)
    println("DCGLP Statistics:")

    dual_eta1 = dual.(env.dcglp.model[:coneta1])
    dual_eta2 = dual.(env.dcglp.model[:coneta2])
    non_zero_indices_dual_eta1 = findall(x -> abs(x) >= 1e-6, dual_eta1)
    non_zero_indices_dual_eta2 = findall(x -> abs(x) >= 1e-6, dual_eta2)
    println("Non-zero dual values for eta1 constraints: ", non_zero_indices_dual_eta1)
    if non_zero_indices_dual_eta1 != []
        println("eta1: ", env.dcglp.model[:coneta1][non_zero_indices_dual_eta1])
        println("value: ", dual_eta1[non_zero_indices_dual_eta1])
    end
    println("Non-zero dual values for eta2 constraints: ", non_zero_indices_dual_eta2)
    if non_zero_indices_dual_eta2 != []
        println("eta2: ", env.dcglp.model[:coneta2][non_zero_indices_dual_eta2])
        println("value: ", dual_eta2[non_zero_indices_dual_eta2])
    end

    dual_conv1 = dual.(env.dcglp.model[:conv1])
    dual_conv2 = dual.(env.dcglp.model[:conv2])
    non_zero_indices_dual_conv1 = findall(x -> abs(x) >= 1e-6, dual_conv1)
    non_zero_indices_dual_conv2 = findall(x -> abs(x) >= 1e-6, dual_conv2)
    println("Non-zero dual values for conv1 constraints: ", non_zero_indices_dual_conv1)
    if non_zero_indices_dual_conv1 != []
        println("conv1: ", env.dcglp.model[:conv1][non_zero_indices_dual_conv1])
        println("value: ", dual_conv1[non_zero_indices_dual_conv1])
    end
    println("Non-zero dual values for conv2 constraints: ", non_zero_indices_dual_conv2)
    if non_zero_indices_dual_conv2 != []
        println("conv2: ", env.dcglp.model[:conv2][non_zero_indices_dual_conv2])
        println("value: ", dual_conv2[non_zero_indices_dual_conv2])
    end

    dual_conk0 = dual.(env.dcglp.model[:conk0])
    dual_conv0 = dual.(env.dcglp.model[:conv0])
    println("Dual values for conk0 constraints: ", dual_conk0)
    println("Dual values for conv0 constraints: ", dual_conv0)

    for i in 1:length(env.dcglp.dcglp_constraints)
        dual_benders_cuts = dual.(env.dcglp.dcglp_constraints[i])
        non_zero_indices_dual_benders_cuts = findall(x -> abs(x) >= 1e-6, dual_benders_cuts)
        println("Non-zero dual values for $(i)th iteration benders cuts: ", non_zero_indices_dual_benders_cuts)
        if non_zero_indices_dual_benders_cuts != []
            println("cuts: ", env.dcglp.dcglp_constraints[i][non_zero_indices_dual_benders_cuts])
            println("value: ", dual_benders_cuts[non_zero_indices_dual_benders_cuts])
        end
    end
    
    σ₁::Float64 = dual.(env.dcglp.disjunctive_inequalities_constraints[1])
    σ₂::Float64 = dual.(env.dcglp.disjunctive_inequalities_constraints[2])
    println("Dual values for disjunctive inequalities constraints: [σ₁: $σ₁, σ₂: $σ₂]")

    γₜ = dual.(env.dcglp.γ_constraints[:γₜ])
    γ₀ = dual.(env.dcglp.γ_constraints[:γ₀])
    γₓ = dual.(env.dcglp.γ_constraints[:γₓ])
    non_zero_indices_dual_γₓ = findall(x -> abs(x) >= 1e-6, γₓ)
    println("γₓ: ", γₓ)
    println("Non-zero dual values for γₓ constraints: ", non_zero_indices_dual_γₓ)
    non_zero_indices_dual_γₜ = findall(x -> abs(x) >= 1e-6, γₜ)
    println("γₜ: ", γₜ)
    println("Non-zero dual values for γₜ constraints: ", non_zero_indices_dual_γₜ)
    disjunctive_cut = @expression(env.master.model, γ₀ + dot(γₓ, env.master.var[:x]) + dot(γₜ, env.master.var[:t]))
    println("disjunctive cut: ", disjunctive_cut)

end

function add_norm_cons(env::BendersEnv)

    # @info "num_constraints of dcglp before deletion", length(all_constraints(env.dcglp.model, include_variable_in_set_constraints = false))
    if haskey(env.dcglp.model, :concone)
        delete.(env.dcglp.model, env.dcglp.model[:concone])
        unregister(env.dcglp.model, :concone)
    end
    # @info "num_constraints of dcglp after deletion", length(all_constraints(env.dcglp.model, include_variable_in_set_constraints = false))

    N = env.data.n_facilities
    M = env.data.n_customers
    dim = 1 + N + M

    # scale s_x

    scaling_value = 10000
    # split_index = env.dcglp.split_index[end]
    # scaling_vector = ones(N)
    # scaling_vector[split_index] = scaling_value * scaling_vector[split_index]
    @constraint(env.dcglp.model, concone, [env.dcglp.model[:τ]; scaling_value .* env.dcglp.model[:sx]; env.dcglp.model[:st]] in MOI.NormInfinityCone(dim))
    # @constraint(env.dcglp.model, concone, [env.dcglp.model[:τ]; scaling_vector .* env.dcglp.model[:sx]; env.dcglp.model[:st]] in MOI.NormInfinityCone(dim))
    # # @constraint(env.dcglp.model, concone, [env.dcglp.model[:τ]; scaling_vector .* env.dcglp.model[:sx]; env.dcglp.model[:st]] in MOI.NormOneCone(dim))

    # # scale s_t

    # scaling_value = 0.00005
    # scaling_vector = scaling_value * ones(N)

    # @constraint(env.dcglp.model, concone, [env.dcglp.model[:τ]; env.dcglp.model[:sx]; scaling_value .* env.dcglp.model[:st]] in MOI.NormInfinityCone(dim))
    # @constraint(env.dcglp.model, concone, [env.dcglp.model[:τ]; env.dcglp.model[:sx]; scaling_value .* env.dcglp.model[:st]] in MOI.NormOneCone(dim))

    # # @info "num_constraints of dcglp after add", length(all_constraints(env.dcglp.model, include_variable_in_set_constraints = false))
end
