export BendersBnBParam

abstract type AbstractBendersBnBParam end

abstract type AbstractBendersBnBState end

abstract type AbstractBendersBnBLog end


mutable struct BendersBnBState <: AbstractBendersBnBState
    oracle_time::Float64
    values::Dict{Symbol,Vector{Float64}}
    f_x::Vector{Float64}
    is_in_L::Bool
    node::Int
    num_cuts::Int
    function BendersBnBState()
        new(0.0, Dict(:x => Vector{Float64}(), :t => Vector{Float64}()), Vector{Float64}(), false, 0, 0)
    end
end

mutable struct BendersBnBLog <: AbstractBendersBnBLog
    nodes::Vector{BendersBnBState}
    n_enter_nodes::Int
    n_lazy_cuts::Int
    n_user_cuts::Int
    start_time::Float64
    num_of_fraction_node::Int
    function BendersBnBLog()
        new(Vector{BendersBnBState}(), 0, 0, 0, 0)
    end
end

"""
Parameters for configuring the Callback-based Benders decomposition algorithm.

Contains settings for:
- `time_limit`: Maximum runtime allowed for the algorithm in seconds.
- `gap_tolerance`: Relative optimality gap tolerance for termination.
- `verbose`: Controls the level of logging output during execution.

These parameters allow fine-tuning of the Benders algorithm performance.
"""
mutable struct BendersBnBParam <: AbstractBendersBnBParam
    time_limit::Float64
    gap_tolerance::Float64
    verbose::Bool

    function BendersBnBParam(; 
                        time_limit::Float64 = 7200.0, 
                        gap_tolerance::Float64 = 1e-6, 
                        verbose::Bool = true
                        ) 
        new(time_limit, gap_tolerance, verbose)
    end
end 

function record_node!(log::BendersBnBLog, state::BendersBnBState, is_lazy_cut::Bool)
    push!(log.nodes, state)
    log.n_enter_nodes += 1
    log.n_lazy_cuts += is_lazy_cut ? state.num_cuts : 0
    log.n_user_cuts += !is_lazy_cut ? state.num_cuts : 0
end