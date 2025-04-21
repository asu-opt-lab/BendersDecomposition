export BendersSeqInOutParam 

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
Check termination criteria for the Sequential Benders loop.

Terminates if:
- `is_in_L` is true (termination via feasibility).
- `gap` is within `gap_tolerance`.
- The remaining time is exhausted.

Returns a `Bool`.
"""
function is_terminated(state::BendersSeqState, log::BendersSeqLog, param::BendersSeqInOutParam)
    return state.is_in_L || state.gap <= param.gap_tolerance || get_sec_remaining(log, param) <= 0.0
end
