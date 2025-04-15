export BendersSeqParam 

abstract type AbstractBendersSeqState <: AbstractLoopState end
abstract type AbstractBendersSeqLog <: AbstractLoopLog end
abstract type AbstractBendersSeqParam <: AbstractLoopParam end
"""
Concrete state for a single iteration of a Sequential Benders algorithm.

Stores:
- `master_time`: Time spent solving the master problem.
- `oracle_time`: Time spent generating cuts by the oracle.
- `total_time`: Combined time for the full iteration.
- `values`: Candidate solution from the master problem.
- `f_x`: Evaluation of the subproblem at `values`; set to `NaN` if the oracle does not evaluate it.
- `is_in_L`: Boolean flag indicating whether the candidate solution is feasible (lies in the feasible region).
- `LB`: Current lower bound.
- `UB`: Current upper bound.
- `gap`: Relative optimality gap (in percent).

Used to store and track relevant data for each iteration of the algorithm.
"""
mutable struct BendersSeqState <: AbstractBendersSeqState
    master_time::Float64
    oracle_time::Float64
    total_time::Float64
    values::Dict{Symbol, Vector{Float64}}
    f_x::Vector{Float64}
    is_in_L::Bool
    LB::Float64
    UB::Float64
    gap::Float64
   
    # Constructor with specified values
    function BendersSeqState()
        new(0.0, 0.0, 0.0, Dict(:x => Vector{Float64}(), :t => Vector{Float64}()), Vector{Float64}(), false, -Inf, Inf, 100.0)
    end
end

"""
Concrete log type for tracking the progress of a Sequential Benders algorithm.

Tracks:
- `n_iter`: Number of iterations completed.
- `iterations`: A vector of `BendersSeqState` objects recording iteration history.
- `start_time`: Timestamp when the algorithm started.
- `consecutive_no_improvement`: Number of iterations with insufficient LB improvement.

Used to store all historical data during the solve and check termination logic.
"""
mutable struct BendersSeqLog <: AbstractBendersSeqLog
    n_iter::Int
    iterations::Vector{BendersSeqState}
    start_time::Float64
    consecutive_no_improvement::Int
    
    function BendersSeqLog()
        new(0, Vector{BendersSeqState}(), time(), 0)
    end
end

mutable struct BendersSeqParam <: AbstractBendersSeqParam

    time_limit::Float64
    gap_tolerance::Float64
    halt_limit::Int
    iter_limit::Int
    verbose::Bool

    function BendersSeqParam(; 
                        time_limit::Float64 = 7200.0, 
                        gap_tolerance::Float64 = 1e-4, 
                        halt_limit::Int = 10000, 
                        iter_limit::Int = 1000000, 
                        verbose::Bool = true
                        ) 
        
        new(time_limit, gap_tolerance, halt_limit, iter_limit, verbose)
    end
end
function update_upper_bound_and_gap!(state::BendersSeqState, log::BendersSeqLog, f::Function)
    evaluation = f(state.f_x, state.values[:x])
    state.UB = log.n_iter >= 1 ? min(log.iterations[end].UB, evaluation) : evaluation
    state.gap = (state.UB - state.LB) / abs(state.UB) * 100
end

function print_iteration_info(state::BendersSeqState, log::BendersSeqLog)
    @printf("Iter: %4d | LB: %12.4f | UB: %11.4f | Gap: %8.3f%% | Time: (M: %6.2f, S: %6.2f, Total: %6.2f) \n",
           log.n_iter, state.LB, state.UB, state.gap, 
           state.master_time, state.oracle_time, state.total_time)
end

"""
Check termination criteria for the Sequential Benders loop.

Terminates if:
- `is_in_L` is true (termination via feasibility).
- `gap` is within `gap_tolerance`.
- The remaining time is exhausted.

Returns a `Bool`.
"""
function is_terminated(state::BendersSeqState, log::BendersSeqLog, param::BendersSeqParam)
    return state.is_in_L || state.gap <= param.gap_tolerance || get_sec_remaining(log, param) <= 0.0
end
