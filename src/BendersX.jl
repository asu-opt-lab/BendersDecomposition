module BendersX

using JuMP
using MathOptInterface
using Printf
using LinearAlgebra
using SparseArrays
using DataFrames

const MOI = MathOptInterface

include("types.jl")
include("utils/utils.jl")
include("modules/modules.jl")
include("artifact_utils.jl")
include("problems/problems.jl")

export AbstractBendersEnv
export AbstractData
export AbstractMaster
export AbstractOracle, AbstractOracleParam, AbstractTypicalOracle, BasicOracleParam
export AbstractBendersSeq, AbstractBendersBnB
export AbstractRootNodePreprocessing, root_node_processing!
export TerminationStatus, NotSolved, TimeLimit, Optimal, InfeasibleOrNumericalIssue
export GBCBoundType, UpperBound, LowerBound, Fixed
export TimeLimitException, UnexpectedModelStatusException, UndefError, AlgorithmException, UnsupportedModelException
export AbstractLoopState, AbstractLoopLog, AbstractLoopParam
export AbstractBendersSeqState, AbstractBendersSeqLog, AbstractBendersSeqParam
export AbstractBendersBnBState, AbstractBendersBnBLog, AbstractBendersBnBParam
export BendersSeqState, BendersSeqLog, BendersSeqParam
export BendersSeqInOutParam, SpecializedBendersSeqParam
export BendersBnBState, BendersBnBLog, BendersBnBParam
export Hyperplane, aggregate, generate_cuts, set_parameter!, hyperplanes_to_expression, select_top_fraction, evaluate_violation, add_constraints
export get_sec_remaining, record_iteration!, update_upper_bound_and_gap!, is_terminated, check_lb_improvement!, print_iteration_info, to_dataframe

export AbstractNorm, StandardNorm, LpNorm
export SplitIndexSelectionRule, RandomFractional, MostFractional, LargestFractional
export DisjunctiveCutsAppendRule, NoDisjunctiveCuts, AllDisjunctiveCuts, DisjunctiveCutsSmallerIndices

end # module BendersX
