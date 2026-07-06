module BendersX

using JuMP
using MathOptInterface
using Printf
using LinearAlgebra
using SparseArrays
using DataFrames
using GLPK

const MOI = MathOptInterface

include("types.jl")
include("utils/utils.jl")
include("modules/modules.jl")
include("artifact_utils.jl")
include("problems/problems.jl")

# Public API surface is centralized here. Child files may define symbols, but
# they must not decide visibility via `export` or `public`.

# Default `using BendersX` workflow
export AbstractData
export Master
export BendersSeq, BendersSeqInOut, SpecializedBendersSeq, BendersBnB
export ClassicalOracle, UnifiedOracle, ParetoOracle, SeparableOracle
export SplitOracle
export solve!

# User-facing parameter and configuration types
export BasicOracleParam, ClassicalOracleParam, UnifiedOracleParam, ParetoOracleParam
export SeparableOracleParam, DcglpParam
export SplitOracleParam
export LpDistanceNormalization
export ReversePolarNormalization
export BendersSeqParam, BendersSeqInOutParam, BendersBnBParam, SpecializedBendersSeqParam
export LazyCallback, UserCallback, NoUserCallback, UserCallbackParam
export RootNodePreprocessing, NoRootNodePreprocessing, DisjunctiveRootNodePreprocessing

# Common runtime statuses and configuration values
export TerminationStatus, NotSolved, TimeLimit, Optimal, InfeasibleOrNumericalIssue
export GBCBoundType, UpperBound, LowerBound, Fixed
export RandomFractional, MostFractional, LargestFractional
export NoDisjunctiveCuts, AllDisjunctiveCuts, DisjunctiveCutsSmallerIndices

# Public extension interfaces
export AbstractBendersEnv, AbstractBendersSeq, AbstractBendersBnB
export AbstractMaster
export AbstractOracle, AbstractOracleParam, AbstractTypicalOracle, AbstractDisjunctiveOracle, AbstractSplitOracle
export AbstractRootNodePreprocessing
export AbstractDisjunctiveNormalization
export SplitIndexSelectionRule, DisjunctiveCutsAppendRule
export generate_cuts, add_normalization_constraint!, update_dcglp_upper_bound_and_gap!
export disjunctive_cut_normalization_value
export callback_node_count, callback_node_depth
export set_parameter!
export set_core_point!
export update_master_model!, update_sub_model!

# Public but not auto-imported advanced utilities and support types
export Hyperplane, aggregate, evaluate_violation, select_top_fraction
export hyperplanes_to_expression, add_constraints
export copy_variables!, var_from_tuple, transfer_scaled_linear_rows_and_bounds_with_types!
export infeasibility_report
export TimeLimitException, UnexpectedModelStatusException, UndefError
export AlgorithmException, UnsupportedModelException

# Problem-specific helpers available via `using BendersX`
export CFLPData, UFLPData, SCFLPData, SNIPData
export CFLKnapsackOracle, CFLKnapsackOracleParam
export UFLKnapsackOracle, UFLKnapsackOracleParam
export read_GK_data, read_cfl_file, read_cflp_benchmark_data
export read_uflp_benchmark_data, read_Simple_data
export read_stochastic_capacited_facility_location_problem, read_snip_data

end # module BendersX
