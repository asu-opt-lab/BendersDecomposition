"""
Utility module for Benders decomposition algorithm implementation.
Provides data structures and helper functions to track algorithm state,
iteration information, and performance metrics.
"""

"""
Stores information about a single iteration of the Benders algorithm.

Fields:
- iter::Int : Current iteration number
- LB::Float64 : Lower bound value
- UB::Float64 : Upper bound value
- gap::Float64 : Optimality gap (%)
- master_time::Float64 : Time spent solving master problem
- sub_time::Float64 : Time spent solving subproblem
- total_time::Float64 : Total elapsed time
"""
struct IterationInfo
    iter::Int
    LB::Float64
    UB::Float64
    gap::Float64
    master_time::Float64
    sub_time::Float64
    total_time::Float64
end

"""
Maintains the current state of the Benders decomposition algorithm.

Fields:
- iteration::Int : Current iteration count
- UB::Float64 : Current upper bound
- LB::Float64 : Current lower bound
- gap::Float64 : Current optimality gap (%)
"""
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

"""
Tracks and stores iteration history and timing information.

Fields:
- iterations::Vector{IterationInfo} : Log of all iterations
- start_time::Float64 : Algorithm start timestamp
- master_time::Float64 : Cumulative master problem solve time
- sub_time::Float64 : Cumulative subproblem solve time
"""
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
    state.UB = min(state.UB, sub_obj_val + dot(env.data.fixed_costs, env.master.x_value))  # UB should always take the minimum value
    update_gap!(state)
end


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


