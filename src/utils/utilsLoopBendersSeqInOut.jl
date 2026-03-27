
"""
    BendersSeqInOutParam(; time_limit = 7200.0, gap_tolerance = 1e-4,
                           halt_limit = 10000, iter_limit = 1000000,
                           verbose = true, stabilizing_x, α = 0.9, λ = 0.1)

Parameter container for [`BendersSeqInOut`](@ref).

In addition to the standard loop controls, this parameter set stores the
initial stabilizing point and the weights used by the in-out stabilization
updates.

# Fields
- `time_limit::Float64`: Global time limit in seconds.
- `gap_tolerance::Float64`: Relative optimality gap tolerance in percent.
- `halt_limit::Int`: Maximum number of consecutive iterations without progress.
- `iter_limit::Int`: Hard iteration limit.
- `verbose::Bool`: Whether to print iteration logs.
- `stabilizing_x::Vector{Float64}`: Initial stabilization point in the master-variable space.
- `α::Float64`: Weight used when updating the stabilization point.
- `λ::Float64`: Weight used when forming the perturbed query point.
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
