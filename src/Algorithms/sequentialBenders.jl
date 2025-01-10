"""
Solve the Benders decomposition algorithm
"""
function solve!(env::BendersEnv, ::Sequential, cut_strategy::CutStrategy, params::BendersParams)
    log = BendersIterationLog()
    state = BendersState()
    
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
            solve_sub!(env.sub, env.master.x_value)
            cuts, sub_obj_val = generate_cuts(env, cut_strategy)
            update_upper_bound_and_gap!(state, env, sub_obj_val)
        end
        log.sub_time += sub_time
        
        # Update state and record information
        record_iteration!(log, state)

        params.verbose && print_iteration_info(state, log)

        # Check termination criteria
        is_terminated(state, params, log) && break

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
    # return state.gap <= params.gap_tolerance || get_total_time(log) >= params.time_limit #|| state.iteration >= 80
    return state.gap <= params.gap_tolerance || get_total_time(log) >= 100 #|| state.iteration >= 80
end