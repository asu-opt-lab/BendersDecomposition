"""
Solve the Benders decomposition algorithm
"""
function solve!(env::BendersEnv, ::Sequential, cut_strategy::CutStrategy, params::BendersParams)
    log = BendersIterationLog()
    state = BendersState()
    stabilizing_x = ones(length(env.master.var[:x]))
    α = 0.9
    λ = 0.1
    relax_integrality(env.master.model)
    prev_lb = -Inf
    consecutive_no_improvement = 0
    kelley_mode = false
    
    while true
        state.iteration += 1
        
        # Solve master problem
        master_time = @elapsed begin
            solve_master!(env.master)
            state.LB = env.master.obj_value
        end
        log.master_time += master_time
        
        # Solve sub problem
        sub_time = @elapsed begin
            stabilizing_x = α * stabilizing_x + (1 - α) * env.master.x_value
            intermediate_x = λ * env.master.x_value + (1 - λ) * stabilizing_x
            solve_sub!(env.sub, intermediate_x)
            cuts, sub_obj_val = generate_cuts(env, cut_strategy, intermediate_x)
            # update_upper_bound_and_gap!(state, env, sub_obj_val)
        end
        log.sub_time += sub_time
        
        # Update state and record information
        record_iteration!(log, state)

        params.verbose && print_iteration_info(state, log)

        # Check termination criteria
        # is_terminated(state, params, log) && break
        
        if state.iteration > 1 && prev_lb != -Inf
            lb_improvement = abs((state.LB - prev_lb) / prev_lb) * 100
            
            # Check for improvement
            if lb_improvement < 0.05
                consecutive_no_improvement += 1
                
                # After 5 consecutive iterations without improvement
                if consecutive_no_improvement == 5 && !kelley_mode
                    # Reset λ to 1 (switch to Kelley's cutting plane)
                    λ = 1.0
                    kelley_mode = true
                    consecutive_no_improvement = 0
                    params.verbose && println("Switching to Kelley's cutting plane method (λ = 1.0)")
                # After 5 more consecutive iterations without improvement in Kelley mode
                elseif consecutive_no_improvement == 5 && kelley_mode
                    params.verbose && println("Aborting after 5 consecutive iterations without improvement in Kelley mode")
                    break
                end
            else
                # Reset counter if there's improvement
                consecutive_no_improvement = 0
            end
        end
        
        prev_lb = state.LB

        # Generate and add cuts
        for cut in cuts
            @constraint(env.master.model, 0 >= cut)
        end
    end
    
    return to_dataframe(log)
end

"""
Solve the master problem
"""
function solve_master!(master::AbstractMasterProblem)
    optimize!(master.model)
    # @info termination_status(master.model)
    master.obj_value = JuMP.objective_value(master.model)
    master.x_value = value.(master.var[:x])
    master.t_value = value.(master.var[:t])
end

"""
Solve the sub problem
"""
function solve_sub!(sub::AbstractSubProblem, x_value::Vector{Float64})
    set_normalized_rhs.(sub.fixed_x_constraints, x_value)
    optimize!(sub.model)
end

"""
Special case for knapsack sub problem - no need to solve
"""
function solve_sub!(sub::KnapsackUFLPSubProblem, x_value::Vector{Float64})
    # No need to solve sub problem (knapsack)
end


"""
Print iteration information if verbose mode is on
"""
function print_iteration_info(state::BendersState, log::BendersIterationLog)
    @printf("Iter: %4d | LB: %12.4f | UB: %11.4f | Gap: %8.3f%% | Time: (M: %6.2f, S: %6.2f) | Elapsed: %6.2f\n",
           state.iteration, state.LB, state.UB, state.gap, 
           log.master_time, log.sub_time, get_total_time(log))
end

"""
Check termination criteria based on gap and time limit
"""
function is_terminated(state::BendersState, params::BendersParams, log::BendersIterationLog)
    return state.gap <= params.gap_tolerance || get_total_time(log) >= params.time_limit #|| state.iteration >= 10
end

function generate_cuts(env::BendersEnv, ::ClassicalCut, intermediate_x::Vector{Float64})
    (coefficients_t, coefficients_x, constant_term), sub_obj_val = generate_cut_coefficients(env.sub, intermediate_x, ClassicalCut())

    cut = @expression(env.master.model, 
        constant_term + dot(coefficients_x, env.master.var[:x]) + coefficients_t * env.master.var[:t])
    # @info "Classical cut: $cut"
    return cut, sub_obj_val
end

function generate_cuts(env::BendersEnv, ::KnapsackCut, intermediate_x::Vector{Float64})
    (μ, KP_values, coeff_t), sub_obj_val = generate_cut_coefficients(env.sub, intermediate_x, KnapsackCut())

    cut = @expression(env.master.model, 
        coeff_t * env.master.var[:t] + μ_term(env.sub, μ) + dot(KP_values, env.master.var[:x]))
    # @info "Knapsack cut: $cut"
    return cut, sub_obj_val
end 