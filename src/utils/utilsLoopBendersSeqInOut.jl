export BendersSeqInOutParam 

"""
    BendersSeqInOutParam <: AbstractBendersSeqParam

Parameter container for the *In–Out* variant of the sequential Benders
decomposition algorithm.

`BendersSeqInOutParam` controls termination criteria, logging behavior, and
stabilization settings used in the sequential In–Out Benders loop. In
particular, it stores a user-provided stabilization point and parameters that
govern how strongly the algorithm is biased toward this reference solution.

# Constructor

```julia
BendersSeqInOutParam(;
    time_limit::Float64 = 7200.0,
    gap_tolerance::Float64 = 1e-4,
    halt_limit::Int = 10000,
    iter_limit::Int = 1_000_000,
    verbose::Bool = true,
    stabilizing_x::Vector{Float64},
    α::Float64 = 0.9,
    λ::Float64 = 0.1,
)
```

# Fields
In addition to the basic fields inherited from `BendersSeqParam` (e.g., `time_limit`, `gap_tolerance`, `halt_limit`, `iter_limit`, `verbose`), this environment introduces the following In–Out–specific fields:
- `stabilizing_x`::Vector{Float64}: User-specified stabilization point for the In–Out algorithm. This vector must be provided and should match the dimension of the master decision variables.
- `α::Float64`: Weight parameter controlling the influence of the stabilizing point (typically 0 < α ≤ 1).
- `λ::Float64`: Step-size or relaxation parameter used in the In–Out update rule.

# Notes
- `stabilizing_x` is mandatory and has no default value.
- The parameters `α` and `λ` jointly control stability versus progress; poor choices may slow convergence.

See also: [`BendersSeqInOut`](@ref), [`BendersSeqParam`](@ref)
"""
mutable struct BendersSeqInOutParam <: AbstractBendersSeqParam

    time_limit::Float64
    gap_tolerance::Float64
    halt_limit::Int
    iter_limit::Int
    verbose::Bool

    stabilizing_x::Vector{Float64} 
    α::Float64
    λ::Float64

    function BendersSeqInOutParam(; 
                        time_limit::Float64 = 7200.0, 
                        gap_tolerance::Float64 = 1e-4, 
                        halt_limit::Int = 10000, 
                        iter_limit::Int = 1000000, 
                        verbose::Bool = true,
                        stabilizing_x::Vector{Float64}, # must be provided
                        α::Float64 = 0.9,
                        λ::Float64 = 0.1
                        ) 
        
        new(time_limit, gap_tolerance, halt_limit, iter_limit, verbose, stabilizing_x, α, λ)
    end
end

"""
Check termination criteria for the Sequential Benders InOut loop.

Terminates if:
- `is_in_L` is true (termination via feasibility).
- `gap` is within `gap_tolerance`.
- The remaining time is exhausted.

Returns a `Bool`.
"""
function is_terminated(state::BendersSeqState, log::BendersSeqLog, param::BendersSeqInOutParam)
    return state.is_in_L || state.gap <= param.gap_tolerance || get_sec_remaining(log, param) <= 0.0
end