"""
Run BendersSeq
"""
function solve!(env::BendersEnv, ::Seq, params::BendersParams) 
try    
    log = BendersIterationLog()
    state = BendersState()
    
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
        
        # Execute oracle
        oracle_time = @elapsed begin
            log.is_in_L, hyperplanes, sub_obj_val = generate_cuts(env.oracle, env.master.x_value, env.master.t_value)

            cuts = !log.is_in_L ? @expression(env.master.model, [j=1:length(hyperplanes)], hyperplanes[j].a_0 + hyperplanes[j].a_x'*env.master.model[:x] + hyperplanes[j].a_t'*env.master.model[:t]) : []

            if sub_obj_val != NaN
                update_upper_bound_and_gap!(state, env, sub_obj_val)
            end
        end
        log.sub_time += oracle_time
        
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
catch e
    @info e
    # 1. generate_cuts throws error when dcglp master is not optimally solved
    # 2. generate_cuts throws error when model-based oracle is not optimally solved
    # 3. it throws error when master is not optimally solved
end
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
    return log.is_in_L || get_total_time(log) >= params.time_limit #|| state.iteration >= 10
    # return state.gap <= params.gap_tolerance || get_total_time(log) >= params.time_limit #|| state.iteration >= 10
end