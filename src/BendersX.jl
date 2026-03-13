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

# Public API surface is centralized here. Child files may define symbols, but
# they must not decide visibility via `export` or `public`.

# Default `using BendersX` workflow
export AbstractData
export Master
export BendersSeq, BendersSeqInOut, SpecializedBendersSeq, BendersBnB
export ClassicalOracle, UnifiedOracle, ParetoOracle, SeparableOracle, SplitOracle
export solve!

# User-facing parameter and configuration types
export BasicOracleParam, ClassicalOracleParam, UnifiedOracleParam, ParetoOracleParam
export SeparableOracleParam, SplitOracleParam, DcglpParam
export BendersSeqParam, BendersSeqInOutParam, BendersBnBParam, SpecializedBendersSeqParam
export LazyCallback, UserCallback, NoUserCallback, UserCallbackParam
export RootNodePreprocessing, NoRootNodePreprocessing, DisjunctiveRootNodePreprocessing

# Common runtime statuses and configuration values
export TerminationStatus, NotSolved, TimeLimit, Optimal, InfeasibleOrNumericalIssue
export GBCBoundType, UpperBound, LowerBound, Fixed
export LpNorm
export RandomFractional, MostFractional, LargestFractional
export NoDisjunctiveCuts, AllDisjunctiveCuts, DisjunctiveCutsSmallerIndices

# Public but not auto-imported extension interfaces
public AbstractBendersEnv, AbstractBendersSeq, AbstractBendersBnB
public AbstractMaster
public AbstractOracle, AbstractOracleParam, AbstractTypicalOracle
public AbstractRootNodePreprocessing
public AbstractNorm, StandardNorm
public SplitIndexSelectionRule, DisjunctiveCutsAppendRule
public generate_cuts, set_parameter!
public customize_master_model!, customize_sub_model!, customize_mip_model!

# Public but not auto-imported advanced utilities and support types
public Hyperplane, aggregate, evaluate_violation, select_top_fraction
public hyperplanes_to_expression, add_constraints
public copy_variables!, var_from_tuple, transfer_scaled_linear_rows_and_bounds_with_types!
public assign_attributes!, infeasibility_report
public TimeLimitException, UnexpectedModelStatusException, UndefError
public AlgorithmException, UnsupportedModelException

# Public but not auto-imported problem-specific helpers
public CFLPData, UFLPData, SCFLPData, SNIPData
public CFLKnapsackOracle, CFLKnapsackOracleParam
public UFLKnapsackOracle, UFLKnapsackOracleParam
public read_GK_data, read_cfl_file, read_cflp_benchmark_data
public read_uflp_benchmark_data, read_Simple_data
public read_stochastic_capacited_facility_location_problem, read_snip_data

end # module BendersX
