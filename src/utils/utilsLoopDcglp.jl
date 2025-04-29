export DcglpParam

abstract type AbstractDcglpState <: AbstractLoopState end
abstract type AbstractDcglpLog <: AbstractLoopLog end
abstract type AbstractDcglpParam <: AbstractLoopParam end

mutable struct DcglpState <: AbstractDcglpState
    master_time::Float64
    oracle_times::Vector{Float64}
    total_time::Float64
    values::Dict{Symbol,Any}
    f_x::Vector{Vector{Float64}}
    omega_t_::Vector{Vector{Float64}}
    is_in_L::Vector{Bool}
    LB::Float64
    UB::Float64
    gap::Float64

    # Constructor with default values
    function DcglpState() 
        new(0.0, 
            [0.0; 0.0], 
            0.0,
            Dict(:ω_x => Vector{Vector{Float64}}(undef, 2), 
                 :ω_t => Vector{Vector{Float64}}(undef, 2), 
                 :ω_0 => Vector{Float64}(undef, 2), 
                 :tau => -Inf,
                 :sx => Vector{Float64}()),
            Vector{Vector{Float64}}(undef, 2), 
            Vector{Vector{Float64}}(undef, 2), 
            [false; false], 
            -Inf, 
            Inf, 
            100.0)
    end
end
"""
Log for DCGLP cutting-plane run
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

mutable struct DcglpParam <: AbstractDcglpParam
    
    time_limit::Float64
    gap_tolerance::Float64
    halt_limit::Int
    iter_limit::Int
    verbose::Bool

    function DcglpParam(; 
                        time_limit::Float64 = 1000.0, 
                        gap_tolerance::Float64 = 1e-3, 
                        halt_limit::Int = 3, 
                        iter_limit::Int = 250,
                        verbose::Bool = true
                        ) 
        
        new(time_limit, gap_tolerance, halt_limit, iter_limit, verbose)
    end
end

"""
utility functions for Dcglp
"""
function update_upper_bound_and_gap!(state::DcglpState, log::DcglpLog, f::Function)
    for i=1:2 
        state.omega_t_[i] = state.is_in_L[i] ? state.values[:ω_t][i] : state.f_x[i] * state.values[:ω_0][i]
    end
    evaluation = f(state.omega_t_[1], state.omega_t_[2])
    
    state.UB = log.n_iter >= 1 ? min(log.iterations[end].UB, evaluation) : evaluation
    state.gap = (state.UB - state.LB) / abs(state.UB) * 100
    # (t1, t2) -> LinearAlgebra.norm([state.other_values[:sx]; t1 .+ t2 .- t_value], norm.p)
end

"""
Print iteration information if verbose mode is on
"""
function print_iteration_info(state::DcglpState, log::DcglpLog)
    @printf("   Iter: %4d | LB: %8.4f | UB: %8.4f | Gap: %6.2f%% | UB_k: %8.2f | UB_v: %8.2f | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
           log.n_iter, state.LB, state.UB, state.gap, sum(state.omega_t_[1]), sum(state.omega_t_[2]), state.master_time, state.oracle_times[1], state.oracle_times[2])
end

"""
Print iteration information if verbose mode is on
"""
function print_iteration_info(state::DcglpState, log::DcglpLog, ::StandardNorm)
    @printf("   Iter: %4d | LB: %8.4f | UB: %8.4f | Gap: %6.2f%% | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
           log.n_iter, state.LB, state.UB, state.gap, state.master_time, state.oracle_times[1], state.oracle_times[2])
end

"""
Check termination criteria
"""
function is_terminated(state::DcglpState, log::DcglpLog, param::DcglpParam, time_limit::Float64)
    return state.is_in_L[1] && state.is_in_L[2] || log.consecutive_no_improvement >= param.halt_limit || state.gap <= param.gap_tolerance  || get_sec_remaining(log.start_time, time_limit) <= 0.0 || time() - log.start_time >= param.time_limit || log.n_iter >= param.iter_limit
end