struct IterationInfo
    iter::Int
    LB::Float64
    UB::Float64
    gap::Float64
    master_time::Float64
    sub_time::Float64
    total_time::Float64
end
# BendersState needed? 
mutable struct BendersState
    iteration::Int
    UB::Float64
    LB::Float64
    gap::Float64

    # Constructor with default values
    function BendersState()
        new(0, Inf, -Inf, Inf)
    end

    # Constructor with specified values
    function BendersState(LB::Float64, UB::Float64, gap::Float64, iteration::Int = 0)
        new(iteration, UB, LB, gap)
    end
end

mutable struct BendersIterationLog
    iterations::Vector{IterationInfo}
    start_time::Float64
    master_time::Float64
    sub_time::Float64
    is_in_L::Bool

    function BendersIterationLog()
        new(IterationInfo[], time(), 0.0, 0.0, false)
    end
end

# Helper functions for BendersIterationLog
function get_total_time(log::BendersIterationLog)
    return time() - log.start_time
end

function record_iteration!(log::BendersIterationLog, state::BendersState)
    push!(log.iterations, IterationInfo(
        state.iteration,
        state.LB,
        state.UB,
        state.gap,
        log.master_time,
        log.sub_time,
        get_total_time(log)
    ))
end

# Helper functions for BendersState
function update_gap!(state::BendersState)
    state.gap = (state.UB - state.LB) / abs(state.UB) * 100
end


function update_upper_bound_and_gap!(state::BendersState, env::BendersEnv, sub_obj_val::Vector{Float64})

    state.UB = min(state.UB, env.data.c_t' * sub_obj_val + env.data.c_x' * env.master.x_value)

    update_gap!(state)
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


function to_dataframe(log::BendersIterationLog)
    return DataFrame(
        iter = [info.iter for info in log.iterations],
        LB = [info.LB for info in log.iterations],
        UB = [info.UB for info in log.iterations],
        gap = [info.gap for info in log.iterations],
        master_time = [info.master_time for info in log.iterations],
        sub_time = [info.sub_time for info in log.iterations],
        total_time = [info.total_time for info in log.iterations]
    )
end

