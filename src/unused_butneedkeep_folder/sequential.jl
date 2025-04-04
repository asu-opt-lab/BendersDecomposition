"""
Solve the Benders decomposition algorithm
"""
function solve!(env::AbstractBendersEnv, SolutionProcedure::GenericSequential, cut_strategy::AbstractCutStrategy)
    # parameters
    time_limit = SolutionProcedure.time_limit
    iteration_limit = SolutionProcedure.iteration_limit
    gap_tolerance = SolutionProcedure.gap_tolerance
    verbose = SolutionProcedure.verbose

    state = BendersState()
    start_time = time()
    while true
        state.iteration += 1
        
        # Solve master problem
        master_time = @elapsed begin
            solve_master!(env.master)
            state.LB = env.master.objective_value
        end
        state.master_time += master_time
        
        # Solve sub problem
        sub_time = @elapsed begin
            solve_sub!(env.sub, env.master.integer_variable_values)
            cuts, sub_obj_val = generate_cuts(env, cut_strategy)
            update_upper_bound_and_gap!(state, env, sub_obj_val)
        end
        state.sub_time += sub_time
        state.total_time = time() - start_time

        # Update state and record information
        verbose && print_iteration_info(state)

        # Check termination criteria
        is_terminated(state, iteration_limit, time_limit, gap_tolerance) && break

        # Generate and add cuts
        add_cuts!(env, cuts)
    end
    
end


function solve!(env::BendersEnv, SolutionProcedure::SequentialWithInout, cut_strategy::CutStrategy)
    # parameters
    time_limit = SolutionProcedure.base.time_limit
    iteration_limit = SolutionProcedure.base.iteration_limit
    gap_tolerance = SolutionProcedure.base.gap_tolerance
    verbose = SolutionProcedure.base.verbose
    stabilizing_point = ones(length(env.master.var[:x]))
    α = SolutionProcedure.α
    λ = SolutionProcedure.λ
    prev_lb = -Inf
    consecutive_no_improvement = 0
    kelley_mode = false
    state = BendersState()
    start_time = time()
    while true
        state.iteration += 1
        
        # Solve master problem
        master_time = @elapsed begin
            solve_master!(env.master)
            state.LB = env.master.objective_value
        end
        state.master_time += master_time
        
        # Solve sub problem
        sub_time = @elapsed begin
            stabilizing_point = α * stabilizing_point + (1 - α) * env.master.integer_variable_values
            intermediate_point = λ * env.master.integer_variable_values + (1 - λ) * stabilizing_point
            solve_sub!(env.sub, intermediate_point)
            cuts, sub_obj_val = generate_cuts(env, cut_strategy, intermediate_point)
            update_upper_bound_and_gap!(state, env, sub_obj_val)
        end
        state.sub_time += sub_time
        state.total_time = time() - start_time

        # Update state and record information
        verbose && print_iteration_info(state)

        is_terminated(state, iteration_limit, time_limit, gap_tolerance) && break

        if prev_lb != -Inf
            lb_improvement = abs((state.LB - prev_lb) / prev_lb) * 100
            
            # Check for improvement
            if lb_improvement < 0.001
                consecutive_no_improvement += 1
                # After 5 consecutive iterations without improvement
                if consecutive_no_improvement == 5 && !kelley_mode
                    # Reset λ to 1 (switch to Kelley's cutting plane)
                    λ = 1.0
                    kelley_mode = true
                    consecutive_no_improvement = 0
                    verbose && println("Switching to Kelley's cutting plane method (λ = 1.0)")
                # After 5 more consecutive iterations without improvement in Kelley mode
                elseif consecutive_no_improvement == 5 && kelley_mode
                    verbose && println("Aborting after 5 consecutive iterations without improvement in Kelley mode")
                    break
                end
            else
                # Reset counter if there's improvement
                consecutive_no_improvement = 0
            end
        end
        prev_lb = state.LB
        # Generate and add cuts
        add_cuts!(env, cuts)
    end
    
end
