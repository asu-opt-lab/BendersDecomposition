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
- `values::Dict{Symbol,Vector{Float64}}`: Solution values at this node, typically includes `:x` (first-stage variables) and `:t` (a vector of auxiliary variables for subproblem objective approximations).
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
- `root_node_time::Float64`: Time spent preprocessing at the root node.
- `total_time::Float64`: Total time spent in `solve!`.
- `num_of_fraction_node::Int`: Number of nodes with fractional solutions.
"""
mutable struct BendersBnBLog <: AbstractBendersBnBLog
    nodes::Vector{BendersBnBState}
    n_enter_nodes::Int
    n_lazy_cuts::Int
    n_user_cuts::Int
    start_time::Float64
    root_node_time::Float64
    total_time::Float64
    num_of_fraction_node::Int
    function BendersBnBLog()
        new(Vector{BendersBnBState}(), 0, 0, 0, 0.0, 0.0, 0.0, 0)
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

"""
Convert the BnB log into a DataFrame.
"""
function to_dataframe(env::AbstractBendersBnB, log::BendersBnBLog)
    if termination_status(env.master.model) == MOI.OPTIMIZE_NOT_CALLED
        return DataFrame(
            node_count = 0,
            root_node_time = log.root_node_time,
            time = log.total_time,
            obj_bound = -Inf,
            obj_val = Inf,
            rel_gap = NaN,
            n_lazy_cuts = log.n_lazy_cuts,
            n_user_cuts = log.n_user_cuts
        )
    else
        return DataFrame(
            node_count = JuMP.node_count(env.master.model),
            root_node_time = log.root_node_time,
            time = log.total_time,
            obj_bound = JuMP.objective_bound(env.master.model),
            obj_val = env.obj_value,
            rel_gap = has_values(env.master.model) ? JuMP.relative_gap(env.master.model) : NaN,
            n_lazy_cuts = log.n_lazy_cuts,
            n_user_cuts = log.n_user_cuts
        )
    end
end