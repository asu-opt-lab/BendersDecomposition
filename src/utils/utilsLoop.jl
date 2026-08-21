"""
    AbstractLoopState

Abstract supertype for state information maintained during an optimization loop.

Concrete subtypes store the current bounds, solution values, subproblem evaluations, and timing information needed to evaluate and record the progress of one iteration.

Loop states are stored in the iteration history of an [`AbstractLoopLog`](@ref).

# Common Fields

Concrete loop-state types typically provide:

- `LB::Float64`: Current lower bound.
- `UB::Float64`: Current upper bound.
- `gap::Float64`: Optimality gap between the upper and lower bounds, in percent.
- `values`: Current solution estimates.
- `f_x`: Subproblem objective values associated with the current solution.
- `is_in_L::Bool`: Whether the current solution belongs to the feasible region of the oracle.
- `master_time`, `oracle_time`, `total_time`: Time spent for each per iteration (optional).

See also: [`AbstractLoopLog`](@ref), [`AbstractLoopParam`](@ref)
"""
abstract type AbstractLoopState end

"""
    AbstractLoopLog

Abstract supertype for the iteration log of an optimization loop.

An `AbstractLoopLog` stores the sequence of [`AbstractLoopState`](@ref) objects generated during the loop, together with aggregate information such as the iteration count, start time, and counters used by termination criteria.

In particular, `iterations` contains the state recorded for each completed iteration of the corresponding loop.

# Common Fields

- `iterations`: History of [`AbstractLoopState`](@ref) objects.
- `n_iter::Int`: Number of completed iterations.
- `start_time::Float64`: Time at which the loop started.
- `consecutive_no_improvement::Int`: Number of consecutive iterations without
  sufficient lower-bound improvement.

See also: [`AbstractLoopState`](@ref), [`AbstractLoopParam`](@ref)
"""

abstract type AbstractLoopLog end

"""
    AbstractLoopParam

Abstract supertype for parameter containers used by optimization loops.

Concrete subtypes provide common loop controls such as time limits, optimality gap tolerances, iteration limits, and verbosity, together with algorithm-specific parameters.

# Common Fields

- `time_limit::Float64`: Maximum wall-clock time allowed for the loop, in
  seconds.
- `gap_tolerance::Float64`: Optimality-gap tolerance used for termination.
- `iter_limit::Int`: Maximum number of loop iterations.
- `verbose::Bool`: Whether to print iteration information.
"""
abstract type AbstractLoopParam end

"""
    print_iteration_info(
        state::AbstractLoopState,
        log::AbstractLoopLog,
    )

Generic interface for printing iteration information.

Concrete loop implementations should specialize this method to report iteration statistics appropriate to their state and log types.

# Throws

Throws [`UnimplementedInterfaceException`](@ref) if no implementation is defined for the supplied state and log types.
"""
function print_iteration_info(
    state::AbstractLoopState,
    log::AbstractLoopLog,
)
    throw(
        UnimplementedInterfaceException(
            "print_iteration_info is not implemented for " *
            "$(typeof(state)) and $(typeof(log)).",
        ),
    )
end


"""
    update_upper_bound_and_gap!(
        state::AbstractLoopState,
        log::AbstractLoopLog,
        f::Function,
    )

Generic interface for updating the upper bound and optimality gap stored in the loop state.

Concrete implementations use `f` to evaluate a candidate upper bound and update `state.UB` and `state.gap`. The gap is reported as a percentage relative to the absolute upper bound.

# Throws

Throws [`UnimplementedInterfaceException`](@ref) if no implementation is defined for the supplied state and log types.
"""
function update_upper_bound_and_gap!(
    state::AbstractLoopState,
    log::AbstractLoopLog,
    f::Function,
)
    throw(
        UnimplementedInterfaceException(
            "update_upper_bound_and_gap! is not implemented for " *
            "$(typeof(state)) and $(typeof(log)).",
        ),
    )
end

"""
    is_terminated(
        state::AbstractLoopState,
        log::AbstractLoopLog,
        param::AbstractLoopParam,
    )

Generic interface for checking whether an optimization loop should terminate.

Concrete implementations define the termination criteria appropriate for the corresponding state, log, and parameter types.

# Throws

Throws [`UnimplementedInterfaceException`](@ref) if no specialized termination rule is defined for the supplied types.
"""
function is_terminated(
    state::AbstractLoopState,
    log::AbstractLoopLog,
    param::AbstractLoopParam,
)
    throw(
        UnimplementedInterfaceException(
            "is_terminated is not implemented for " *
            "$(typeof(state)), $(typeof(log)), and $(typeof(param)).",
        ),
    )
end

"""
    check_lb_improvement!(
        state::AbstractLoopState,
        log::AbstractLoopLog;
        zero_tol = 1e-8,
        tol_imprv = 1e-4,
    )

Check the improvement in the lower bound and update the consecutive no-improvement counter.

The improvement is measured as the absolute change in the lower bound when the previous lower bound is within `zero_tol` of zero, and as the relative change otherwise.

If the measured improvement is smaller than `tol_imprv`, `log.consecutive_no_improvement` is incremented. Otherwise, the counter is reset to zero.

# Arguments

- `state::AbstractLoopState`: Current loop state containing the lower bound.
- `log::AbstractLoopLog`: Loop log containing previous iterations and the consecutive no-improvement counter.
- `zero_tol`: Tolerance below which the previous lower bound is treated as zero.
- `tol_imprv`: Minimum lower-bound improvement required to reset the no-improvement counter.
"""
function check_lb_improvement!(
    state::AbstractLoopState,
    log::AbstractLoopLog;
    zero_tol = 1e-8,
    tol_imprv = 1e-4,
)
    prev_lb = log.n_iter > 1 ? log.iterations[end - 1].LB : state.LB

    lb_improvement =
        abs(prev_lb) < zero_tol ?
        abs(state.LB - prev_lb) :
        abs((state.LB - prev_lb) / prev_lb)

    if lb_improvement < tol_imprv
        log.consecutive_no_improvement += 1
    else
        log.consecutive_no_improvement = 0
    end

    return nothing
end

"""
    get_sec_remaining(start_time::Float64, time_limit::Float64)

Return the number of seconds remaining in a wall-clock time budget.

The returned value is clamped to zero when the time limit has been reached.

See also: [`AbstractLoopLog`](@ref)
"""
function get_sec_remaining(tic::Float64, time_limit::Float64)
    time_elapsed = time() - tic
    return max(time_limit - time_elapsed, 0.0)
end

"""
    get_sec_remaining(log::AbstractLoopLog, param::AbstractLoopParam)

Return the number of seconds remaining in the loop's configured time budget.

The remaining time is computed from `log.start_time` and `param.time_limit`.

See also: [`get_sec_remaining`](@ref)
"""
function get_sec_remaining(log::AbstractLoopLog, param::AbstractLoopParam)
    return get_sec_remaining(log.start_time, param.time_limit)
end

"""
    record_iteration!(
        log::AbstractLoopLog,
        state::AbstractLoopState,
    )

Append the current iteration state to the loop log and increment the iteration counter.

Returns `nothing`.
"""
function record_iteration!(log::AbstractLoopLog, state::AbstractLoopState)
    push!(log.iterations, state)
    log.n_iter += 1
end

"""
    to_dataframe(log::AbstractLoopLog) -> DataFrame

Convert the iteration log into a `DataFrame`.

The resulting data frame contains the lower bound, upper bound, optimality gap, and available timing statistics for each iteration.
"""
function to_dataframe(log::AbstractLoopLog)
    return DataFrame(
        LB = [info.LB for info in log.iterations],
        UB = [info.UB for info in log.iterations],
        gap = [info.gap for info in log.iterations],
        master_time = [info.master_time for info in log.iterations],
        oracle_time = [info.oracle_time for info in log.iterations],
        total_time = [info.total_time for info in log.iterations]
    )
end

# ---------------------------------------------------------------------------- 
# Concrete implementations 
# ----------------------------------------------------------------------------
include("utilsLoopBendersSeq.jl")
include("utilsLoopBendersSeqInOut.jl")
include("utilsLoopDcglp.jl")
