export get_sec_remaining

"""
Abstract type representing the state of an algorithm at a given iteration.

Concrete subtypes of `AbstractLoopState` should store the current values of:
- `LB::Float64`: Current lower bound.
- `UB::Float64`: Current upper bound.
- `gap::Float64`: Optimality gap between UB and LB (in percent).
- `values`: Current solution estimate.
- `f_x`: Subproblem objective value at `values`; set to `NaN` if the oracle does not evaluate it.
- `master_time`, `oracle_time`, `total_time`: Timing statistics per iteration (optional).

This type is used to track progress and update algorithmic decisions at each iteration.
"""
abstract type AbstractLoopState end

"""
Abstract type representing the log of iterations during the optimization loop.

Concrete subtypes should store:
- `iterations::Vector{<:AbstractLoopState}`: History of all iteration states.
- `n_iter::Int`: Total number of iterations performed.
- `start_time::Float64`: Time when the algorithm started (used for time limits).
- `consecutive_no_improvement::Int`: Counter for how many iterations have passed without LB improvement.

This type is useful for performance tracking, debugging, and termination logic.
"""
abstract type AbstractLoopLog end

"""
Abstract type for parameter containers used in the optimization loop.

Concrete subtypes should store:
- `time_limit::Float64`: Time budget for the algorithm in seconds.
- `gap_tolerance::Float64`: .
- `iter_limit::Float64`: .
- `verbose::Bool`: Whether and how much to log to stdout.
- Other algorithm-specific parameters for controlling behavior and tuning.
"""
abstract type AbstractLoopParam end

"""
Print iteration information if verbose mode is on.
This is a prototype and can be implemented for specific subtypes of `AbstractLoopState` and `AbstractLoopLog`.
"""
function print_iteration_info(state::AbstractLoopState, log::AbstractLoopLog)
    throw(UndefError("update print_iteration_info for $(typeof(state)) and $(typeof(log))"))
end

"""
Update the upper bound and optimality gap in the loop state: state.UB, state.gap. The UB is updated based on a user-supplied evaluation function `f`. The gap is computed as a percentage relative to the absolute UB.
This is a prototype and can be implemented for specific subtypes of `AbstractLoopState` and `AbstractLoopLog`.
"""
function update_upper_bound_and_gap!(state::AbstractLoopState, log::AbstractLoopLog, f::Function)
    throw(UndefError("update update_upper_bound_and_gap! for $(typeof(state)) and $(typeof(log))"))
end

"""
Check whether the loop should terminate.
This is a prototype and can be implemented for specific subtypes of `AbstractLoopState` and `AbstractLoopLog`.
"""
function is_terminated(state::AbstractLoopState, log::AbstractLoopLog, params::AbstractLoopParam)
    throw(UndefError("update is_terminated for $(typeof(state)) and $(typeof(log))"))
end

"""
Check for improvement in the lower bound (LB) and update the no-improvement counter.
If the LB improvement is smaller than `tol_imprv` (default 1e-4), increment the `consecutive_no_improvement` counter. Otherwise, reset the counter.
"""
function check_lb_improvement!(state::AbstractLoopState, log::AbstractLoopLog; zero_tol = 1e-8, tol_imprv = 1e-4)
    prev_lb = log.n_iter > 1 ? log.iterations[end-1].LB : state.LB
    lb_improvement = abs(prev_lb) < zero_tol ? abs(state.LB - prev_lb) : abs((state.LB - prev_lb) / prev_lb) * 100
    # Check for improvement
    if lb_improvement < tol_imprv
        log.consecutive_no_improvement += 1
    else
        # Reset counter if there's improvement
        log.consecutive_no_improvement = 0
    end
end

"""
Append constraints to a JuMP model using a symbolic name.
If constraints with the name `constr_symbol` exist, append to them. Otherwise, create a new constraint group named `constr_symbol`, enforcing `0 .>= exprs`.
"""
function add_constraints(model::Model, constr_symbol::Symbol, exprs::Vector{AffExpr})
    # add constraints in the form of 0 .>= expr
    if haskey(model, constr_symbol)
        append!(model[constr_symbol], @constraint(model, 0 .>= exprs))
    else
        model[constr_symbol] = @constraint(model, 0 .>= exprs)
    end
end

"""
Return the number of seconds remaining given a start time and time limit.
Ensures the returned time is no smaller than `tol` (default 1e-4).
"""
function get_sec_remaining(tic::Float64, time_limit::Float64; tol = 1e-4)
    # tol needed to prevent parameter being too small
    time_elapsed = time() - tic
    return max(time_limit - time_elapsed, tol)
end
"""
Return the number of seconds remaining using the start time given in log and the time limit specified in a parameter object.
"""
function get_sec_remaining(log::AbstractLoopLog, param::AbstractLoopParam)
    return get_sec_remaining(log.start_time, param.time_limit)
end

"""
Record the current iteration state into the loop log.
"""
function record_iteration!(log::AbstractLoopLog, state::AbstractLoopState)
    push!(log.iterations, state)
    log.n_iter += 1
end

"""
Convert the loop log into a DataFrame.
"""
function to_dataframe(log::AbstractLoopLog)
    return DataFrame(
        LB = [info.LB for info in log.iterations],
        UB = [info.UB for info in log.iterations],
        gap = [info.gap for info in log.iterations],
        master_time = [info.master_time for info in log.iterations],
        oracle_time = [info.oracle_time for info in log.iterations],
        total_time = [info.total_time for info in log.iterations],
        # is_in_L = [info.is_in_L for info in log.iterations]
    )
end

include("utilsLoopBendersSeq.jl")
include("utilsLoopDcglp.jl")
include("utilsLoopBendersSeqInOut.jl")
include("utilsLoopSpecializedBendersSeq.jl")