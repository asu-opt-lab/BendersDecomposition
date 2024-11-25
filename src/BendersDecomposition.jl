module BendersDecomposition

using Printf, StatsBase, Random, Distributions, LinearAlgebra, Plots, ArgParse, DataFrames, CSV, JSON
using JuMP, CPLEX, Gurobi

# Include supporting files
include("types.jl")
include("models/models.jl")
include("utils/utils.jl")


"""
    BendersParams

Configuration parameters for the Benders Decomposition algorithm.

# Fields
- `time_limit::Float64`: Maximum execution time
- `gap_tolerance::Float64`: Convergence tolerance
- `main_verbose::Bool`: Enable main algorithm logging
- `dcglp_verbose::Bool`: Enable DCGLP logging
- `solver::Symbol`: Selected optimization solver
- `master_attributes::Dict{Symbol,Any}`: Master problem solver attributes
- `sub_attributes::Dict{Symbol,Any}`: Sub problem solver attributes
- `dcglp_attributes::Dict{Symbol,Any}`: DCGLP solver attributes
"""
struct BendersParams
    time_limit::Float64
    gap_tolerance::Float64
    solver::Symbol
    master_attributes::Dict{Symbol,Any}
    sub_attributes::Dict{Symbol,Any}
    dcglp_attributes::Dict{Symbol,Any}
    verbose::Bool
end
export BendersParams

"""
    BendersEnv

Main structure for the Benders Decomposition algorithm.

# Fields
- `data::D`: Problem data
- `loop_strategy::L`: Strategy for main algorithm loop
- `cut_strategy::C`: Strategy for generating Benders cuts
- `master::AbstractMasterProblem`: Master problem formulation
- `sub::AbstractSubProblem`: Sub problem formulation
- `dcglp::Union{Nothing, DCGLP}`: Optional DCGLP component
"""
mutable struct BendersEnv
    data::AbstractData
    master::AbstractMasterProblem
    sub::AbstractSubProblem
    dcglp::Union{Nothing, DCGLP}  # Optional component
end

function BendersEnv(data::AbstractData, cut_strategy::CutStrategy
, params::BendersParams)
    master = create_master_problem(data, cut_strategy)
    assign_attributes!(master.model, params.master_attributes)
    sub = create_sub_problem(data, cut_strategy)
    assign_attributes!(sub.model, params.sub_attributes)
    if cut_strategy isa DisjunctiveCut
        dcglp = create_dcglp(data, cut_strategy)
        assign_attributes!(dcglp.model, params.dcglp_attributes)
    else
        dcglp = nothing
    end
    return BendersEnv(data, master, sub, dcglp)
end
export BendersEnv

"""
    run_Benders(data::D, loop_strategy::L, cut_strategy::C, params::BendersParams)

Execute the Benders Decomposition algorithm with the given configuration.

# Arguments
- `data::D`: Problem instance data
- `loop_strategy::L`: Strategy for main algorithm loop
- `cut_strategy::C`: Strategy for generating Benders cuts
- `params::BendersParams`: Algorithm parameters

# Returns
- Solution from the Benders Decomposition algorithm
"""
function run_Benders(data::AbstractData, loop_strategy::SolutionProcedure, cut_strategy::CutStrategy, params::BendersParams)
    env = BendersEnv(data, cut_strategy, params)
    solve!(env, loop_strategy, cut_strategy, params)
end
export run_Benders

include("algorithms/algorithms.jl")

end
