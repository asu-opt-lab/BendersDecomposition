struct IterationInfo
    iter::Int
    LB::Float64
    UB::Float64
    gap::Float64
    master_time::Float64
    sub_time::Float64
    total_time::Float64
end

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

    function BendersIterationLog()
        new(IterationInfo[], time(), 0.0, 0.0)
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


function update_upper_bound_and_gap!(state::BendersState, env::BendersEnv, sub_obj_val::Float64)

    state.UB = min(state.UB, sub_obj_val + dot(env.data.fixed_costs, env.master.x_value))

    update_gap!(state)
end


function update_upper_bound_and_gap!(state::BendersState, env::BendersEnv, sub_obj_val_collection::Vector{Float64})

    state.UB = min(state.UB, sum(sub_obj_val_collection[k] * env.data.scenarios[k][3] for k in 1:env.data.num_scenarios))  
    update_gap!(state)
end




