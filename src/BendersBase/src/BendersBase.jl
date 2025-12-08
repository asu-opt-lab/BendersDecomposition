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
export Data, AbstractData, EmptyData
export AbstractMaster, AbstractMip
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
# 1. Remove mip.jl and Mip struct
# 2. Remove EmptyData; remove AbstractData completely and consider accepting Any
# 3. remove previous initializer of master, oracle, env modules
# 4. remove Data by having dim_x, dim_t, c_x, c_t info in master module
# 4. rename oracle_param to param for all oracles as we no longer has solver_param