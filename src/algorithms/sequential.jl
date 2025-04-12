"""
Run BendersSeq
"""
function solve!(env::BendersEnv, ::Seq, params::BendersParams) 
log = BendersLog()
try    
    state = BendersState()
    while true
        state.iteration += 1
        
        # Solve master problem
        state.master_time = @elapsed begin
            set_time_limit_sec(env.master.model, get_sec_remaining(log, params))
            optimize!(env.master.model)
            if is_solved_and_feasible(env.master.model; allow_local = false, dual = false)
                env.master.obj_value = JuMP.objective_value(env.master.model)
                env.master.x_value = value.(env.master.model[:x])
                env.master.t_value = value.(env.master.model[:t])
                state.LB = env.master.obj_value
            elseif termination_status(env.master.model) == TIME_LIMIT
                throw(TimeLimitException("Time limit reached during master solving"))
            else 
                throw(UnexpectedModelStatusException("Master: $(termination_status(env.master.model))"))
                # if infeasible, then the milp is infeasible
            end
        end
        # log.master_time += state.master_time
        
        # Execute oracle
        state.oracle_time = @elapsed begin
            state.is_in_L, hyperplanes, sub_obj_val = generate_cuts(env.oracle, env.master.x_value, env.master.t_value; time_limit = get_sec_remaining(log, params))
            
            cuts = !state.is_in_L ? hyperplanes_to_expression(env.master.model, hyperplanes, env.master.model[:x], env.master.model[:t]) : []
        
            if sub_obj_val !== NaN
                update_upper_bound_and_gap!(state, env, sub_obj_val)
            end
        end
        # log.oracle_time += state.oracle_time
        
        # Update state and record information
        record_iteration!(log, state)

        params.verbose && print_iteration_info(state, log)

        # Check termination criteria
        is_terminated(state, params, log) && break

        # Generate and add cuts
        @constraint(env.master.model, 0 .>= cuts)
    end
    env.log.termination_status = Optimal()
    
    return to_dataframe(log)
catch e
    if typeof(e) <: TimeLimitException
        log.termination_status = TimeLimit()
    elseif typeof(e) <: UnexpectedModelStatusException
        log.termination_status = InfeasibleOrNumericalIssue()
    else
        rethrow()  
    end
    env.log = log
end

# even if it terminates in the middle due to time limit, should be able to access the latest x_value via env.master
end

"""
Print iteration information if verbose mode is on
"""
function print_iteration_info(state::BendersState, log::BendersLog)
    @printf("Iter: %4d | LB: %12.4f | UB: %11.4f | Gap: %8.3f%% | Time: (M: %6.2f, S: %6.2f) | Elapsed: %6.2f\n",
           state.iteration, state.LB, state.UB, state.gap, 
           state.master_time, state.oracle_time, time() - log.start_time)
end

"""
Check termination criteria based on gap and time limit
"""
function is_terminated(state::BendersState, params::BendersParams, log::BendersLog)
    return state.is_in_L || state.gap <= params.gap_tolerance || get_sec_remaining(log, params) <= 0.0
end
