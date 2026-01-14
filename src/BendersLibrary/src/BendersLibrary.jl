module BendersLibrary

using JuMP
using BendersBase
using Printf
using SparseArrays

import BendersBase: solve!, generate_cuts, update_upper_bound_and_gap!, is_terminated, print_iteration_info, set_parameter!, customize_master_model!, customize_sub_model!, root_node_processing!

include("types.jl")
include("utils/utils.jl")
include("modules/modules.jl") 
include("artifact_utils.jl")
include("problems/problems.jl")

export AbstractNorm, StandardNorm, LpNorm
export SplitIndexSelectionRule, RandomFractional, MostFractional, LargestFractional
export DisjunctiveCutsAppendRule, NoDisjunctiveCuts, AllDisjunctiveCuts, DisjunctiveCutsSmallerIndices

# Re-export solver utilities from BendersBase
export get_cplex_optimizer, get_cplex_lp_optimizer
export get_gurobi_optimizer, get_gurobi_lp_optimizer
export get_available_solvers, is_solver_available

end