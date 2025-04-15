
abstract type AbstractBendersSeq <: AbstractBendersDecomposition end
abstract type AbstractBendersCallback <: AbstractBendersDecomposition end

include("BendersSeq.jl") 
include("BendersSeqInOut.jl") 
include("Dcglp.jl") 

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
