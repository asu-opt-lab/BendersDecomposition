"""
    solve!(env::BendersEnv, loop::SolutionProcedure, cut_strategy::CutStrategy, params::BendersParams)

Main entry point for solving a problem using Benders decomposition.

# Arguments
- `env::BendersEnv`: Environment containing problem data and solution state
- `loop::SolutionProcedure`: Solution procedure type (Sequential or Callback)
- `cut_strategy::CutStrategy`: Strategy for generating Benders cuts
- `params::BendersParams`: Algorithm parameters and settings

# Returns
- Solution status and optimal values

# Throws
- `MethodError`: If the concrete type does not implement required methods
- `ArgumentError`: If input parameters are invalid
"""
# function solve!(env::BendersEnv, loop::SolutionProcedure, cut_strategy::CutStrategy, params::BendersParams)
#     error("solve! not implemented for $(typeof(loop)) with $(typeof(cut_strategy))")
# end

"""
    generate_cuts(env::BendersEnv, cut_strategy::CutStrategy)

Generate Benders cuts based on the current master solution.

# Arguments
- `env::BendersEnv`: Environment containing current solution state
- `cut_strategy::CutStrategy`: Strategy for generating cuts

# Returns
- Collection of generated Benders cuts

# Throws
- `ArgumentError`: If the strategy type is not supported
"""
# function generate_cuts(env::BendersEnv, cut_strategy::CutStrategy)
#     error("generate_cuts not implemented for strategy type $(typeof(cut_strategy))")
# end

"""
    generate_cut_coefficients(sub::AbstractSubProblem, x_value::Vector{Float64}, cut_strategy::CutStrategy)

Generate coefficients for Benders cuts based on subproblem solution.
Must be implemented by concrete cut strategy types.

# Arguments
- `sub::AbstractSubProblem`: Subproblem instance
- `x_value::Vector{Float64}`: Subproblem solution
- `cut_strategy::CutStrategy`: Strategy for generating cuts

# Returns
- Coefficients for Benders cuts

# Throws
- `ArgumentError`: If the subproblem type or strategy type is not supported
"""
# function generate_cut_coefficients(sub::AbstractSubProblem, x_value::Vector{Float64}, cut_strategy::CutStrategy)
#     error("generate_cut_coefficients not implemented for subproblem type $(typeof(sub)) and strategy $(typeof(cut_strategy))")
# end

# Include algorithm implementations
include("utils.jl") 
include("sequential.jl") 
include("sequentialInOut.jl") 

# include("algorithms_utils.jl")
# include("sequentialBenders.jl") 
# include("sequentialBenders_stochastic.jl")
# include("sequentialBenders_test.jl") 
# include("callbackBenders.jl")
# include("callbackBenders_stochastic.jl")

# Include cut strategy implementations
# include("cut_strategies/classical_cut.jl")
# include("cut_strategies/knapsack_cut.jl")
# include("cut_strategies/fs_knapsack_cut.jl")

# # Include disjunctive cut system
# include("disjunction_system/dcglp_cut.jl")
# include("disjunction_system/dcglp_cut_stochastic.jl")
# include("disjunction_system/dcglp_cut_lifting.jl")
# include("disjunction_system/update_dcglp.jl")
# include("disjunction_system/solve_dcglp.jl")
# include("disjunction_system/dcglp_utils.jl")
