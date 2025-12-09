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