abstract type AbstractDcglpState <: AbstractLoopState end
abstract type AbstractDcglpLog <: AbstractLoopLog end
abstract type AbstractDcglpParam <: AbstractLoopParam end

"""
    DcglpState <: AbstractDcglpState

State recorded for one iteration of the DCGLP cutting-plane loop.

# Fields

- `LB::Float64`: Current lower bound.
- `UB::Float64`: Current upper bound.
- `gap::Float64`: Current optimality gap, in percent.
- `values::Dict{Symbol,Any}`: Current DCGLP solution values.
- `f_x::Vector{Vector{Float64}}`: Subproblem objective values obtained from the typical oracles on the two sides of the split.
- `master_time::Float64`: Time spent solving the DCGLP master problem.
- `oracle_times::Vector{Float64}`: Time spent evaluating the typical oracle on each side of the split.
- `total_time::Float64`: Total time spent in the iteration.
- `omega_t_::Vector{Vector{Float64}}`: Auxiliary-vector estimates used to evaluate the current upper bound.
- `is_in_L::Vector{Bool}`: Whether the two disjunctive points belong to the feasible region of their respective typical oracles.

See also: [`DcglpLog`](@ref), [`DcglpParam`](@ref), [`SplitOracle`](@ref)
"""
mutable struct DcglpState <: AbstractDcglpState
    LB::Float64
    UB::Float64
    gap::Float64
    values::Dict{Symbol,Any}
    f_x::Vector{Vector{Float64}}
    is_in_L::Vector{Bool}
    master_time::Float64
    oracle_times::Vector{Float64}
    total_time::Float64
    omega_t_::Vector{Vector{Float64}}
    
    # Constructor with default values
    function DcglpState() 
        new(-Inf, 
            Inf, 
            100.0,
            Dict(:ω_x => Vector{Vector{Float64}}(undef, 2), 
                 :ω_t => Vector{Vector{Float64}}(undef, 2), 
                 :ω_0 => Vector{Float64}(undef, 2), 
                 :tau => -Inf,
                 :sx => Vector{Float64}()),
            Vector{Vector{Float64}}(undef, 2), 
            [false; false], 
            0.0, 
            [0.0; 0.0], 
            0.0,
            Vector{Vector{Float64}}(undef, 2)
            )
    end
end

"""
    DcglpLog <: AbstractDcglpLog

Log of iterations performed by the DCGLP cutting-plane loop.

# Fields

- `n_iter::Int`: Number of completed DCGLP iterations.
- `iterations::Vector{DcglpState}`: Iteration history.
- `start_time::Float64`: Time at which the DCGLP loop started.
- `consecutive_no_improvement::Int`: Number of consecutive iterations without sufficient lower-bound improvement.

See also: [`DcglpState`](@ref), [`DcglpParam`](@ref)
"""
mutable struct DcglpLog <: AbstractDcglpLog
    n_iter::Int
    iterations::Vector{DcglpState}
    start_time::Float64
    consecutive_no_improvement::Int
    
    function DcglpLog()
        new(0, Vector{DcglpState}(), time(), 0)
    end
end

"""
    DcglpParam <: AbstractDcglpParam

Parameters controlling the DCGLP cutting-plane solution process used by [`SplitOracle`](@ref).

# Fields

- `optimizer::MOI.OptimizerWithAttributes`: Optimizer used to solve the DCGLP.
- `time_limit::Float64`: Maximum time allowed for the DCGLP solution process, in seconds.
- `gap_tolerance::Float64`: Optimality-gap tolerance for terminating the DCGLP cutting-plane loop.
- `halt_limit::Int`: Maximum number of consecutive iterations without sufficient lower-bound improvement before termination.
- `iter_limit::Int`: Maximum number of iterations of the DCGLP cutting-plane loop.
- `verbose::Bool`: Whether to print DCGLP iteration information.

The DCGLP cutting-plane loop terminates when a configured stopping criterion is satisfied, including the time limit, gap tolerance, halt limit, or iteration limit.

# Constructor

    DcglpParam(;
        optimizer = default_optimizer(),
        time_limit = 1000.0,
        gap_tolerance = 1e-3,
        halt_limit = 3,
        iter_limit = 250,
        verbose = true,
    )

Construct DCGLP parameters with the specified optimizer and stopping
criteria.

See also: [`SplitOracleParam`](@ref), [`SplitOracle`](@ref)
"""
mutable struct DcglpParam <: AbstractDcglpParam
    
    optimizer::MOI.OptimizerWithAttributes
    time_limit::Float64
    gap_tolerance::Float64
    halt_limit::Int
    iter_limit::Int
    verbose::Bool

    function DcglpParam(;
        optimizer::MOI.OptimizerWithAttributes = default_optimizer(),
        time_limit::Float64 = 1000.0,
        gap_tolerance::Float64 = 1e-3,
        halt_limit::Int = 3,
        iter_limit::Int = 250,
        verbose::Bool = true,
    )
        new(optimizer, time_limit, gap_tolerance, halt_limit, iter_limit, verbose)
    end
end

"""
    update_upper_bound_and_gap!(
        state::DcglpState,
        log::DcglpLog,
        f::Function;
        zero_tol = 1e-9,
    )

Update the upper bound and optimality gap for a DCGLP iteration.

For each side of the split, the auxiliary-vector estimate is taken from the DCGLP solution when the corresponding point is in the oracle feasible region. Otherwise, the estimate is obtained from the subproblem objective values scaled by the corresponding `omega_0` value. 

The function `f` is evaluated at the resulting auxiliary vectors to obtain a candidate upper bound. The best upper bound seen so far is stored in `state.UB`.

The gap is reported as a percentage relative to `abs(state.UB)`. If the upper and lower bounds agree within `zero_tol`, the gap is set to zero. If the upper bound is within `zero_tol` of zero but does not agree with the lower bound, the gap is set to `Inf`.

Returns `nothing`.
"""
function update_upper_bound_and_gap!(
    state::DcglpState,
    log::DcglpLog,
    f::Function;
    zero_tol::Float64 = 1e-9,
)
    for i in 1:2
        state.omega_t_[i] =
            state.is_in_L[i] ?
            state.values[:ω_t][i] :
            state.f_x[i] * state.values[:ω_0][i]
    end

    evaluation = f(
        state.omega_t_[1],
        state.omega_t_[2],
    )

    state.UB =
        log.n_iter >= 1 ?
        min(log.iterations[end].UB, evaluation) :
        evaluation

    if isapprox(state.UB, state.LB; atol = zero_tol)
        state.gap = 0.0
    elseif abs(state.UB) <= zero_tol
        state.gap = Inf
    else
        state.gap = max(
            0.0,
            (state.UB - state.LB) / abs(state.UB) * 100.0,
        )
    end

    return nothing
end

"""
    print_iteration_info(
        state::DcglpState,
        log::DcglpLog,
    )

Print iteration statistics for the DCGLP cutting-plane loop.

The output includes the lower bound, upper bound, optimality gap, auxiliary upper-bound estimates, and master/oracle timing information.
"""
function print_iteration_info(state::DcglpState, log::DcglpLog)
    @printf("   Iter: %4d | LB: %8.4f | UB: %8.4f | Gap: %6.2f%% | UB_k: %8.2f | UB_v: %8.2f | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
           log.n_iter, state.LB, state.UB, state.gap, sum(state.omega_t_[1]), sum(state.omega_t_[2]), state.master_time, state.oracle_times[1], state.oracle_times[2])
end

"""
    is_terminated(
        state::DcglpState,
        log::DcglpLog,
        param::DcglpParam
    )

Check whether the DCGLP cutting-plane loop should terminate.

The loop terminates when any of the following conditions is satisfied:

- both disjunctive branches are in the oracle feasible region;
- the current optimality gap is within `param.gap_tolerance`;
- the lower bound has failed to improve for `param.halt_limit` consecutive iterations;
- the number of iterations reaches `param.iter_limit`.

Time-limit termination is handled by the surrounding solve loop.

Returns `true` if the loop should terminate and `false` otherwise.
"""
function is_terminated(state::DcglpState, log::DcglpLog, param::DcglpParam)
    return (
        all(state.is_in_L) ||
        log.consecutive_no_improvement >= param.halt_limit ||
        state.gap <= param.gap_tolerance ||
        log.n_iter >= param.iter_limit
    )
end