module BendersDecomposition

using Printf, StatsBase, Random, Distributions, LinearAlgebra, ArgParse, DataFrames, CSV, JSON
using JuMP, CPLEX #, Gurobi

# Include supporting files
include("types.jl")
include("models/models.jl") 
include("utils/utils.jl")

"""
    BendersParams(time_limit, gap_tolerance, solver, master_attributes, sub_attributes, dcglp_attributes, verbose)

Parameters for configuring the Benders decomposition algorithm.

# Arguments
- `time_limit::Float64`: Maximum allowed runtime in seconds
- `gap_tolerance::Float64`: Optimality gap tolerance for convergence
- `solver::String`: Optimization solver to use (e.g. "CPLEX", "Gurobi")
- `master_attributes::Dict{String,Any}`: Solver-specific attributes for master problem
- `sub_attributes::Dict{String,Any}`: Solver-specific attributes for subproblem
- `dcglp_attributes::Dict{String,Any}`: Solver-specific attributes for DCGLP
- `verbose::Bool`: Whether to print detailed progress information
"""
mutable struct BendersParams
    time_limit::Float64
    gap_tolerance::Float64
    solver::String
    master_attributes::Dict{String,Any}
    sub_attributes::Dict{String,Any}
    dcglp_attributes::Dict{String,Any}
    verbose::Bool
end
export BendersParams

"""
    BendersEnv(data, master, sub, dcglp)

Environment containing all components needed for Benders decomposition.

# Arguments
- `data::AbstractData`: Problem instance data
- `master::AbstractMasterProblem`: Master problem formulation
- `sub::Union{AbstractSubProblem,Vector{AbstractSubProblem}}`: Subproblem formulation(s)
- `dcglp::Union{Nothing,DCGLP}`: Optional DCGLP component for cut generation
"""
mutable struct BendersEnv
    data::AbstractData
    master::AbstractMasterProblem
    sub::Union{AbstractSubProblem, Vector{AbstractSubProblem}}
    dcglp::Union{Nothing, DCGLP}  # Optional component
end

"""
    BendersEnv(data, cut_strategy, params)

Construct a BendersEnv with the given problem data and configuration.

# Arguments
- `data::AbstractData`: Problem instance data
- `cut_strategy::CutStrategy`: Strategy for generating Benders cuts
- `params::BendersParams`: Algorithm parameters
"""
function BendersEnv(data::AbstractData, cut_strategy::CutStrategy, params::BendersParams)
    master = create_and_configure_master(data, cut_strategy, params)
    sub = create_and_configure_sub(data, cut_strategy, params)
    dcglp = create_and_configure_dcglp(data, cut_strategy, params)
    return BendersEnv(data, master, sub, dcglp)
end

# Helper functions
function create_and_configure_master(data, cut_strategy, params)
    master = create_master_problem(data, cut_strategy)
    assign_attributes!(master.model, params.master_attributes)
    return master
end

function create_and_configure_sub(data, cut_strategy, params)
    sub = create_sub_problem(data, cut_strategy)
    if sub isa Union{AbstractSCFLPSubProblem, AbstractSNIPSubProblem}
        foreach(scenario_sub -> assign_attributes!(scenario_sub.model, params.sub_attributes), 
                sub.sub_problems)
    elseif !isa(sub, KnapsackUFLPSubProblem)
        assign_attributes!(sub.model, params.sub_attributes)
    end
    return sub
end

function create_and_configure_dcglp(data, cut_strategy, params)
    cut_strategy isa DisjunctiveCut || return nothing
    dcglp = create_dcglp(data, cut_strategy)
    assign_attributes!(dcglp.model, params.dcglp_attributes)
    return dcglp
end

export BendersEnv

"""
    run_Benders(data, loop_strategy, cut_strategy, params)

Execute Benders decomposition algorithm to solve the given problem instance.

# Arguments
- `data::AbstractData`: Problem instance data
- `loop_strategy::SolutionProcedure`: Strategy for main algorithm loop (Sequential or Callback)
- `cut_strategy::CutStrategy`: Strategy for generating Benders cuts
- `params::BendersParams`: Algorithm parameters

# Returns
- `DataFrame`: Solution statistics including bounds and timing information
"""
function run_Benders(data::AbstractData, loop_strategy::SolutionProcedure, cut_strategy::CutStrategy, params::BendersParams)
    env = BendersEnv(data, cut_strategy, params)
    solve!(env, loop_strategy, cut_strategy, params)
end

export run_Benders

include("algorithms/algorithms.jl")

end
