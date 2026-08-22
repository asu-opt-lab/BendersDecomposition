"""
    AbstractBendersSeqState <: AbstractLoopState

Abstract supertype for iteration state used by sequential Benders environments.

See also: [`BendersSeqState`](@ref)
"""
abstract type AbstractBendersSeqState <: AbstractLoopState end


"""
    AbstractBendersSeqLog <: AbstractLoopLog

Abstract supertype for iteration logs used by sequential Benders environments.

See also: [`BendersSeqLog`](@ref)
"""
abstract type AbstractBendersSeqLog <: AbstractLoopLog end


"""
    AbstractBendersSeqParam <: AbstractLoopParam

Abstract supertype for parameter containers used by sequential Benders environments.

See also: [`BendersSeqParam`](@ref),
[`BendersSeqInOutParam`](@ref)
"""
abstract type AbstractBendersSeqParam <: AbstractLoopParam end

"""
    BendersSeqState <: AbstractBendersSeqState

State recorded for one iteration of a sequential Benders algorithm.

# Fields

- `LB::Float64`: Current lower bound.
- `UB::Float64`: Current upper bound.
- `gap::Float64`: Current relative optimality gap.
- `values::Dict{Symbol,Vector{Float64}}`: Candidate master solution values, including `:x` and `:t`.
- `f_x::Vector{Float64}`: Subproblem objective values associated with the candidate solution. Entries may be `Inf` when the corresponding objective value is unavailable.
- `master_time::Float64`: Time spent solving the master problem.
- `oracle_time::Float64`: Time spent evaluating the oracle and generating cuts.
- `total_time::Float64`: Total time spent in the iteration.
- `is_in_L::Bool`: Whether the candidate solution belongs to the oracle's feasible region.

See also: [`BendersSeqLog`](@ref), [`BendersSeq`](@ref)
"""
mutable struct BendersSeqState <: AbstractBendersSeqState
    LB::Float64
    UB::Float64
    gap::Float64
    values::Dict{Symbol, Vector{Float64}}
    f_x::Vector{Float64}
    is_in_L::Bool
    master_time::Float64
    oracle_time::Float64
    total_time::Float64
    
    # Constructor with specified values
    function BendersSeqState()
        new(-Inf, Inf, 100.0, Dict(:x => Vector{Float64}(), :t => Vector{Float64}()), Vector{Float64}(), false, 0.0, 0.0, 0.0)
    end
end

"""
    BendersSeqLog <: AbstractBendersSeqLog

Log of iterations performed by a sequential Benders environment.

# Fields

- `n_iter::Int`: Number of completed iterations.
- `iterations::Vector{BendersSeqState}`: Iteration history.
- `start_time::Float64`: Time at which the solve started.
- `consecutive_no_improvement::Int`: Number of consecutive iterations without sufficient lower-bound improvement.
- `preprocessing_time::Float64`: Time spent in preprocessing before the main Benders iterations.

See also: [`BendersSeqState`](@ref), [`BendersSeq`](@ref)
"""
mutable struct BendersSeqLog <: AbstractBendersSeqLog
    n_iter::Int
    iterations::Vector{BendersSeqState}
    start_time::Float64
    consecutive_no_improvement::Int
    preprocessing_time::Float64
    
    function BendersSeqLog()
        new(0, Vector{BendersSeqState}(), time(), 0, 0.0)
    end
end

"""
    BendersSeqParam <: AbstractBendersSeqParam

Parameter container for the standard sequential Benders environment.

`BendersSeqParam` controls global stopping criteria and iteration logging for
[`BendersSeq`](@ref).

# Constructor
```julia
BendersSeqParam(;
    time_limit::Float64 = 7200.0,
    gap_tolerance::Float64 = 1e-4,
    halt_limit::Int = 10000,
    iter_limit::Int = 1_000_000,
    verbose::Bool = true,
)
```

# Fields
- `time_limit::Float64`: Maximum wall-clock time, in seconds.
- `gap_tolerance::Float64`: Relative gap tolerance for termination.
- `halt_limit::Int`: Maximum number of consecutive non-improving iterations.
- `iter_limit::Int`: Maximum number of Benders iterations.
- `verbose::Bool`: Whether to print iteration-level logging information.

See also: [`BendersSeq`](@ref)
"""
mutable struct BendersSeqParam <: AbstractBendersSeqParam

    time_limit::Float64
    gap_tolerance::Float64
    halt_limit::Int
    iter_limit::Int
    verbose::Bool

    function BendersSeqParam(; 
                        time_limit::Float64 = 7200.0, 
                        gap_tolerance::Float64 = 1e-6, 
                        halt_limit::Int = 10000, 
                        iter_limit::Int = 1000000, 
                        verbose::Bool = true
                        ) 
        
        new(time_limit, gap_tolerance, halt_limit, iter_limit, verbose)
    end
end

"""
    update_upper_bound_and_gap!(
        state::BendersSeqState,
        log::BendersSeqLog,
        f::Function;
        zero_tol = 1e-9,
    )

Update the upper bound and optimality gap for a sequential Benders iteration.

The candidate upper bound is obtained by evaluating `f` at the componentwise maximum of the current auxiliary-variable value `t` and the corresponding subproblem objective values `f_x`. The best upper bound seen so far is stored in `state.UB`.

The gap is computed as a relative to `abs(state.UB)`. If the upper bound and lower bound agree within `zero_tol`, the gap is set to zero. If the upper bound is within `zero_tol` of zero but does not agree with the lower bound, the gap is set to `Inf`.

Returns `nothing`.
"""
function update_upper_bound_and_gap!(
    state::BendersSeqState,
    log::BendersSeqLog,
    f::Function;
    zero_tol::Float64 = 1e-9,
)
    evaluation = f(
        max.(state.values[:t], state.f_x),
        state.values[:x],
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
            (state.UB - state.LB) / abs(state.UB),
        )
    end

    return nothing
end

function print_iteration_info(state::BendersSeqState, log::BendersSeqLog; prefix="")
    @printf("%s Iter: %4d | LB: %12.4f | UB: %11.4f | Gap: %8.3f%% | Time: (M: %6.2f, S: %6.2f, Total: %6.2f) \n",
           prefix, log.n_iter, state.LB, state.UB, state.gap * 100, 
           state.master_time, state.oracle_time, state.total_time)
end

"""
    is_terminated(
        state::BendersSeqState,
        log::BendersSeqLog,
        param::BendersSeqParam,
    )

Check whether the sequential Benders loop should terminate.

The loop terminates when the current candidate solution is in the oracle feasible region or when the current optimality gap is within `param.gap_tolerance`.

Time-limit handling is performed by the surrounding solve loop.

Returns `true` if the loop should terminate and `false` otherwise.
"""
function is_terminated(state::BendersSeqState, log::BendersSeqLog, param::BendersSeqParam)
    return state.is_in_L || state.gap <= param.gap_tolerance
end
