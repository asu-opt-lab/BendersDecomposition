"""
Solve the Benders decomposition algorithm
"""
function solve!(env::BendersEnv, ::StochasticSequential, cut_strategy::CutStrategy, params::BendersParams)
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
        sub_obj_val_collection = Float64[]
        cuts_collection = []
        sub_time = @elapsed begin
            for scenario in 1:env.data.num_scenarios
                solve_sub!(env.sub.sub_problems[scenario], env.master.x_value)
                cuts, sub_obj_val = generate_cuts(env, cut_strategy, scenario)
                push!(sub_obj_val_collection, sub_obj_val)
                push!(cuts_collection, cuts)
            end
        end
        update_upper_bound_and_gap!(state, env, sub_obj_val_collection)
        log.sub_time += sub_time
        
        # Update state and record information
        record_iteration!(log, state)

        params.verbose && print_iteration_info(state, log)

        # Check termination criteria
        is_terminated(state, params, log) && break

        # Generate and add cuts
        for cuts in cuts_collection
            for cut in cuts
                @constraint(env.master.model, 0 >= cut)
            end
        end
    end
    
    return to_dataframe(log)
end

function solve!(env::BendersEnv, ::StochasticSequential, cut_strategy::DisjunctiveCut, params::BendersParams)
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
            cuts, sub_obj_val = generate_cuts_stochastic(env, cut_strategy)
        end
        # @info sub_obj_val
        update_upper_bound_and_gap!(state, env, sub_obj_val)
        log.sub_time += sub_time
        
        # Update state and record information
        record_iteration!(log, state)

        params.verbose && print_iteration_info(state, log)

        # Check termination criteria
        is_terminated(state, params, log) && break

        # Generate and add cuts
        for cut in cuts
            # @info "cut: $cut"
            @constraint(env.master.model, 0 .>= cut)
        end
        # println(env.master.model)
    end
    
    return to_dataframe(log)
end