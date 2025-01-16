function generate_cuts_stochastic(env::BendersEnv, cut_strategy::DisjunctiveCut)

    sub_obj_val = get_subproblems_value_stochastic(env) #checked

    disjunctive_inequality = select_disjunctive_inequality(env.master.x_value)

    update_dcglp!(env.dcglp, disjunctive_inequality, cut_strategy)
    
    solve_dcglp_stochastic!(env, cut_strategy)
    
    cuts = merge_cuts_stochastic(env, cut_strategy)

    @info "through"
    return cuts, sub_obj_val
end

function get_subproblems_value_stochastic(env::BendersEnv)
    # # Check if master.x_value is close enough to integer values
    if !all(x -> isapprox(x, round(x), atol=1e-4), env.master.x_value)
        return Inf*ones(env.data.num_scenarios)
    end
    sub_obj_val_collection = Float64[]
    for scenario in 1:env.data.num_scenarios
        solve_sub!(env.sub.sub_problems[scenario], env.master.x_value)
        if dual_status(env.sub.sub_problems[scenario].model) == FEASIBLE_POINT
            push!(sub_obj_val_collection, objective_value(env.sub.sub_problems[scenario].model))
        elseif dual_status(env.sub.sub_problems[scenario].model) == INFEASIBLE_POINT
            push!(sub_obj_val_collection, Inf)
        else
            error("Subproblem is not feasible or optimal")
        end
    end
    return sub_obj_val_collection
end

function solve_dcglp_stochastic!(env::BendersEnv, cut_strategy::DisjunctiveCut)

    log = DCGLPIterationLog()
    state = DCGLPState()

    x_value, t_value = env.master.x_value, env.master.t_value
    set_normalized_rhs.(env.dcglp.model[:conx], x_value)
    set_normalized_rhs.(env.dcglp.model[:cont], t_value)
    log.start_time = time()

    while true

        state.iteration += 1
        if_add_cuts_k_collection = []
        if_add_cuts_v_collection = []
        dual_info_k_collection = []
        dual_info_v_collection = []
        obj_value_k_collection = []
        obj_value_v_collection = []
        master_time = @elapsed begin
            k_values, v_values, other_values = solve_and_get_dcglp_values(env.dcglp.model, cut_strategy.norm_type)
        end
        log.master_time += master_time

        sub_k_time = @elapsed begin
            for scenario in 1:env.data.num_scenarios

                if isapprox(k_values.constant, 0.0, atol=1e-05)
                    if_add_cuts_k, dual_info_k, obj_value_k = false, [], k_values.t[scenario]
                else
                    if_add_cuts_k, dual_info_k, obj_value_k = generate_cut_coefficients_stochastic(env.sub.sub_problems[scenario], k_values, cut_strategy.base_cut_strategy)
                    if obj_value_k <= k_values.t[scenario] + 1e-04
                        if_add_cuts_k, dual_info_k, obj_value_k = false, [], k_values.t[scenario]
                    end
                end
                push!(if_add_cuts_k_collection, if_add_cuts_k)
                push!(dual_info_k_collection, dual_info_k)
                push!(obj_value_k_collection, obj_value_k)
                if isapprox(v_values.constant, 0.0, atol=1e-05)
                    if_add_cuts_v, dual_info_v, obj_value_v = false, [], v_values.t[scenario]
                else
                    if_add_cuts_v, dual_info_v, obj_value_v = generate_cut_coefficients_stochastic(env.sub.sub_problems[scenario], v_values, cut_strategy.base_cut_strategy)
                    if obj_value_v <= v_values.t[scenario] + 1e-04
                        if_add_cuts_v, dual_info_v, obj_value_v = false, [], v_values.t[scenario]
                    end
                end
                push!(if_add_cuts_v_collection, if_add_cuts_v)
                push!(dual_info_v_collection, dual_info_v)
                push!(obj_value_v_collection, obj_value_v)
            end
        end
        log.sub_k_time += sub_k_time / 2
        log.sub_v_time += sub_k_time / 2

        update_bounds!(state, k_values, v_values, other_values, obj_value_k_collection, obj_value_v_collection, t_value, cut_strategy.norm_type)
        
        record_iteration!(log, state)

        cut_strategy.verbose && print_dcglp_iteration_info_stochastic(state, log)

        is_terminated_stochastic(state, log) && break
        
        for scenario in 1:env.data.num_scenarios

            if if_add_cuts_k_collection[scenario]
                add_cuts_k!(env, dual_info_k_collection[scenario], cut_strategy, scenario)
            end
            if if_add_cuts_v_collection[scenario]
                add_cuts_v!(env, dual_info_v_collection[scenario], cut_strategy, scenario)
            end
        end
        
    end

end

function generate_cut_coefficients_stochastic(sub::AbstractSubProblem, korv_values, base_cut_strategy::CutStrategy)
    input = @. abs(korv_values.x / korv_values.constant) # abs value
    solve_sub!(sub, input)
    dual_values, obj_value = generate_cut_coefficients(sub, input, base_cut_strategy)
    obj_value *= korv_values.constant
    return true, dual_values, obj_value
end

function print_dcglp_iteration_info_stochastic(state, log)
    @printf("   Iter: %4d | LB: %12.4f | UB: %11.4f | Gap: %8.2f%% \n",
           state.iteration, state.LB, state.UB, state.gap)
end

function is_terminated_stochastic(state, log)
    return state.gap <= 1e-3  || state.UB - state.LB <= 1e-03 || get_total_time(log) >= 200 || state.iteration >= 20
end

# ============================================================================
# Different types of cuts for the DCGLP
# ============================================================================

function add_cuts_k!(env::BendersEnv, dual_info_k, cut_strategy::CutStrategy, scenario::Int)
    cuts_k, cuts_v, cuts_master = build_cuts(env.dcglp, env.master, env.sub.sub_problems[scenario], dual_info_k, cut_strategy.base_cut_strategy, scenario)
    push!(env.dcglp.dcglp_constraints, @constraint(env.dcglp.model, 0 .>= cuts_k))
    if cut_strategy.use_two_sided_cuts 
        push!(env.dcglp.dcglp_constraints, @constraint(env.dcglp.model, 0 .>= cuts_v))
    end
    if cut_strategy.include_master_cuts 
        push!(env.dcglp.master_cuts, cuts_master)
    end
end

function add_cuts_v!(env::BendersEnv, dual_info_v, cut_strategy::CutStrategy, scenario::Int)
    cuts_k, cuts_v, cuts_master = build_cuts(env.dcglp, env.master, env.sub.sub_problems[scenario], dual_info_v, cut_strategy.base_cut_strategy, scenario)
    push!(env.dcglp.dcglp_constraints, @constraint(env.dcglp.model, 0 .>= cuts_v))
    if cut_strategy.use_two_sided_cuts 
        push!(env.dcglp.dcglp_constraints, @constraint(env.dcglp.model, 0 .>= cuts_k))
    end
    if cut_strategy.include_master_cuts 
        push!(env.dcglp.master_cuts, cuts_master)
    end
end

# knapsack cuts

function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, coeff_info::Tuple{Vector{Float64}, Vector{Float64}, Float64}, ::KnapsackCut, scenario::Int)
    μ, KP_values, coeff_t = coeff_info
    cut_k = @expression(dcglp.model, coeff_t * dcglp.model[:kₜ][scenario] + sum(μ) * dcglp.model[:k₀] + sum(KP_values .* dcglp.model[:kₓ]))
    cut_v = @expression(dcglp.model, coeff_t * dcglp.model[:vₜ][scenario] + sum(μ) * dcglp.model[:v₀] + sum(KP_values .* dcglp.model[:vₓ]))
    cut = @expression(master.model, coeff_t * master.var[:t][scenario] + sum(μ) + sum(KP_values .* master.var[:x]))
    return [cut_k], [cut_v], [cut]
end

function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, coeff_info::Tuple{Float64, Vector{Float64}, Float64}, ::KnapsackCut, scenario::Int)
    coefficients_t, coefficients_x, constant_term = coeff_info
    cut_k = @expression(dcglp.model, constant_term*dcglp.model[:k₀] + sum(coefficients_x .* dcglp.model[:kₓ]) + coefficients_t * dcglp.model[:kₜ][scenario])
    cut_v = @expression(dcglp.model, constant_term*dcglp.model[:v₀] + sum(coefficients_x .* dcglp.model[:vₓ]) + coefficients_t * dcglp.model[:vₜ][scenario])
    cut = @expression(master.model, constant_term + sum(coefficients_x .* master.model[:x]) + coefficients_t * master.model[:t][scenario])
    return [cut_k], [cut_v], [cut]
end


# standard cut
function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, coeff_info::Tuple{Float64, Vector{Float64}, Float64}, ::ClassicalCut, scenario::Int)
    coefficients_t, coefficients_x, constant_term = coeff_info
    cut_k = @expression(dcglp.model, constant_term*dcglp.model[:k₀] + sum(coefficients_x .* dcglp.model[:kₓ]) + coefficients_t * dcglp.model[:kₜ][scenario])
    cut_v = @expression(dcglp.model, constant_term*dcglp.model[:v₀] + sum(coefficients_x .* dcglp.model[:vₓ]) + coefficients_t * dcglp.model[:vₜ][scenario])
    cut = @expression(master.model, constant_term + sum(coefficients_x .* master.model[:x]) + coefficients_t * master.model[:t][scenario])
    return [cut_k], [cut_v], [cut]
end

function merge_cuts_stochastic(env::BendersEnv, cut_strategy::DisjunctiveCut)
    γ₀, γₓ, γₜ = generate_disjunctive_cuts(env.dcglp, cut_strategy.cut_strengthening_type)
    push!(env.dcglp.γ_values, (γ₀, γₓ, γₜ))
    master_disjunctive_cut = @expression(env.master.model, γ₀ + dot(γₓ, env.master.var[:x]) + dot(γₜ, env.master.var[:t]))
    if cut_strategy.include_master_cuts
        push!(env.dcglp.master_cuts, [master_disjunctive_cut])
        return env.dcglp.master_cuts
    end
    return master_disjunctive_cut
end