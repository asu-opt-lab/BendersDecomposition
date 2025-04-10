export get_sec_remaining

# Helper functions for BendersDecompositionLog
function get_sec_remaining(tic::Float64, time_limit::Float64; tol = 1e-4)
    # tol needed to prevent parameter being too small
    time_elapsed = time() - tic
    return max(time_limit - time_elapsed, tol)
end
function get_sec_remaining(log::BendersDecompositionLog, param::BendersParams)
    return get_sec_remaining(log.start_time, param.time_limit)
end
function record_iteration!(log::BendersDecompositionLog, state::BendersState)
    push!(log.iterations, state)
    log.LB = state.LB
    log.UB = state.UB
end

# # Helper functions for BendersState
# function update_gap!(state::BendersState)
#     state.gap = (state.UB - state.LB) / abs(state.UB) * 100
# end


function update_upper_bound_and_gap!(state::BendersState, env::BendersEnv, sub_obj_val::Vector{Float64})

    state.UB = min(state.UB, env.data.c_t' * sub_obj_val + env.data.c_x' * env.master.x_value)
    state.gap = (state.UB - state.LB) / abs(state.UB) * 100
end


# function update_upper_bound_and_gap!(state::BendersState, env::BendersEnv, sub_obj_val_collection::Vector{Float64})
#     if env.data isa SCFLPData
#         state.UB = min(state.UB, mean(sub_obj_val_collection) + dot(env.data.c, env.master.x_value))  
#     elseif env.data isa UFLPData
#         state.UB = min(state.UB, sum(sub_obj_val_collection) + dot(env.data.c, env.master.x_value))  
#     else
#         state.UB = min(state.UB, sum(sub_obj_val_collection[k] * env.data.scenarios[k][3] for k in 1:env.data.num_scenarios))  
#     end
#     update_gap!(state)
# end


function to_dataframe(log::BendersDecompositionLog)
    return DataFrame(
        iter = [info.iteration for info in log.iterations],
        LB = [info.LB for info in log.iterations],
        UB = [info.UB for info in log.iterations],
        gap = [info.gap for info in log.iterations],
        master_time = [info.master_time for info in log.iterations],
        oracle_time = [info.oracle_time for info in log.iterations],
        total_time = [info.total_time for info in log.iterations],
        is_in_L = [info.is_in_L for info in log.iterations]
    )
end

