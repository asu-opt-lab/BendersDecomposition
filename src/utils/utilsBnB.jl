export BendersBnBParam

"""
Abstract type for Benders Branch-and-Bound parameters.
Serves as a parent type for concrete parameter implementations.
"""
abstract type AbstractBendersBnBParam end

"""
Abstract type for Benders Branch-and-Bound state.
Defines the interface for tracking the algorithm's state during execution.
"""
abstract type AbstractBendersBnBState end

"""
Abstract type for Benders Branch-and-Bound logging.
Defines the interface for recording algorithm performance and statistics.
"""
abstract type AbstractBendersBnBLog end


"""
Represents the state of a node in the Branch-and-Bound tree for Benders decomposition.

Fields:
- `oracle_time::Float64`: Time spent in oracle calls for this node.
- `values::Dict{Symbol,Vector{Float64}}`: Solution values at this node, typically includes `:x` (first-stage variables) and `:t` (a vector of subproblem objective approximations).
- `f_x::Vector{Float64}`: A vector of subproblem objective function values at this node.
- `is_in_L::Bool`: Flag indicating if the node is in the L-set.
- `node::Int`: Unique identifier for this node.
- `num_cuts::Int`: Number of Benders cuts generated at this node.
"""
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

"""
Logs execution statistics for the Benders Branch-and-Bound algorithm.

Fields:
- `nodes::Vector{BendersBnBState}`: Collection of all node states explored.
- `n_enter_nodes::Int`: Total number of nodes processed.
- `n_lazy_cuts::Int`: Total number of lazy Benders cuts generated.
- `n_user_cuts::Int`: Total number of user Benders cuts generated.
- `start_time::Float64`: Algorithm start time (used for calculating elapsed time).
- `num_of_fraction_node::Int`: Number of nodes with fractional solutions.
"""
mutable struct BendersBnBLog <: AbstractBendersBnBLog
    nodes::Vector{BendersBnBState}
    n_enter_nodes::Int
    n_lazy_cuts::Int
    n_user_cuts::Int
    start_time::Float64
    num_of_fraction_node::Int
    function BendersBnBLog()
        new(Vector{BendersBnBState}(), 0, 0, 0, 0, 0)
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
    disjunctive_root_process::Bool
    no_good_cut::Bool
    verbose::Bool

    function BendersBnBParam(; 
                        time_limit::Float64 = 7200.0, 
                        gap_tolerance::Float64 = 1e-6, 
                        disjunctive_root_process::Bool = false,
                        no_good_cut::Bool = false,
                        verbose::Bool = true
                        ) 
        new(time_limit, gap_tolerance, disjunctive_root_process, no_good_cut, verbose)
    end
end 

"""
Records a node's state in the Branch-and-Bound log.

# Arguments
- `log::BendersBnBLog`: The log to update.
- `state::BendersBnBState`: The node state to record.
- `is_lazy_cut::Bool`: Whether the cuts generated at this node are lazy cuts.

Updates counters for nodes processed and cuts generated (both lazy and user cuts).
"""
function record_node!(log::BendersBnBLog, state::BendersBnBState, is_lazy_cut::Bool)
    push!(log.nodes, state)
    log.n_enter_nodes += 1
    log.n_lazy_cuts += is_lazy_cut ? state.num_cuts : 0
    log.n_user_cuts += !is_lazy_cut ? state.num_cuts : 0
end

function get_sec_remaining(log::BendersBnBLog, param::BendersBnBParam)
    return get_sec_remaining(log.start_time, param.time_limit)
end