module BendersBase

# Import dependencies
using JuMP
using MathOptInterface
using Printf
using LinearAlgebra
using SparseArrays
const MOI = MathOptInterface
using DataFrames

# Include source files
include("types.jl")
include("utils/utils.jl")
include("modules/modules.jl") 
include("algorithms/algorithms.jl")


# Export abstract types
export AbstractBendersDecomposition
export AbstractData
export AbstractMaster
export AbstractOracle, AbstractOracleParam, AbstractTypicalOracle, BasicOracleParam
export AbstractBendersSeq, AbstractBendersCallback
export Seq, SeqInOut
export TerminationStatus, NotSolved, TimeLimit, Optimal, InfeasibleOrNumericalIssue
export TimeLimitException, UnexpectedModelStatusException, UndefError, AlgorithmException
export AbstractLoopState, AbstractLoopLog, AbstractLoopParam
export AbstractBendersSeqState, AbstractBendersSeqLog, AbstractBendersSeqParam
export AbstractBendersBnBState, AbstractBendersBnBLog, AbstractBendersBnBParam
export BendersSeqState, BendersSeqLog, BendersSeqParam
export BendersBnBState, BendersBnBLog, BendersBnBParam
export Hyperplane, aggregate, generate_cuts, set_parameter!, hyperplanes_to_expression, select_top_fraction, evaluate_violation, add_constraints
export get_sec_remaining, record_iteration!, update_upper_bound_and_gap!, is_terminated, check_lb_improvement!, print_iteration_info, to_dataframe

end # module BendersBase

# To-Do: 
# 1. Remove mip.jl and Mip struct --> done
# 2. remove Data by having dim_x, dim_t, c_x, c_t info in master module --> done
# 3. clean initializers of master, oracle, env modules --> done
# 4. rename `problem::AbstractData` to `data::AbstractData`
# 5. rename oracle_param to param for all oracles as we no longer has solver_param
# 6. consider returning `to_dataframe(log)` for all solve! functions.
# 7. Remove AbstractCallbackParam and EmptyCallbackParam <: AbstractCallbackParam and add AbstractUserCallbackParam; Lazy callback does not need parameters.
# 8. Kaiwen: refactor all files in script folder
# 9. Inho: refactor all snip-related files