"""
Run BendersSeqInOut
"""
function solve!(env::BendersEnv, ::SeqInOut, params::BendersParams)
try    
    log = BendersIterationLog()
    state = BendersState()
    stabilizing_x = ones(env.data.dim_x)
    α = 0.9
    λ = 0.1
    # relax_integrality(env.master.model)
    prev_lb = -Inf
    consecutive_no_improvement = 0
    kelley_mode = false
    
    while true
        state.iteration += 1
        
        # Solve master problem
        master_time = @elapsed begin
            optimize!(env.master.model)
            if is_solved_and_feasible(env.master.model; allow_local = false, dual = false)
                env.master.obj_value = JuMP.objective_value(env.master.model)
                env.master.x_value = value.(env.master.model[:x])
                env.master.t_value = value.(env.master.model[:t])
                state.LB = env.master.obj_value
            else 
                throw(ErrorException("master termination status: $(termination_status(env.master.model))"))
                # if infeasible, then the milp is infeasible
            end
        end
        log.master_time += master_time
        
        # perturb point
        stabilizing_x = α * stabilizing_x + (1 - α) * env.master.x_value
        intermediate_x = λ * env.master.x_value + (1 - λ) * stabilizing_x

        # Execute oracle
        oracle_time = @elapsed begin
            
            log.is_in_L, hyperplanes, sub_obj_val = generate_cuts(env.oracle, env.master.x_value, env.master.t_value)

            # need to differential dim_t = 1 or > 1
            cuts = !log.is_in_L ? @expression(env.master.model, [j=1:length(hyperplanes)],
            hyperplanes[j].a_0 + hyperplanes[j].a_x'*env.master.model[:x] + hyperplanes[j].a_t'*env.master.model[:t]) : []

            if kelley_mode && sub_obj_val != NaN
                # Check termination criteria
                update_upper_bound_and_gap!(state, env, sub_obj_val)
            else
                log.is_in_L = false
            end
        end
        log.sub_time += oracle_time

        # Update state and record information
        record_iteration!(log, state)

        params.verbose && print_iteration_info(state, log)

        is_terminated(state, params, log) && break
        
        if !kelley_mode && state.iteration > 1 && prev_lb != -Inf
            lb_improvement = abs(prev_lb) < 1e-4 ? abs(state.LB - prev_lb) : abs((state.LB - prev_lb) / prev_lb) * 100
            
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
                # elseif consecutive_no_improvement == 5 && kelley_mode
                #     params.verbose && println("Aborting after 5 consecutive iterations without improvement in Kelley mode")
                #     break
                end
            else
                # Reset counter if there's improvement
                consecutive_no_improvement = 0
            end
        end
        
        prev_lb = state.LB

        # Generate and add cuts
        @constraint(env.master.model, 0 .>= cuts)
    end
    
    return to_dataframe(log)
catch e
    @info e
    # 1. generate_cuts throws error when dcglp master is not optimally solved
    # 2. generate_cuts throws error when model-based oracle is not optimally solved
    # 3. it throws error when master is not optimally solved
end
end