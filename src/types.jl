# ============================================================================
# Exports
# ============================================================================

export AbstractBendersDecomposition
export Data, AbstractData
export AbstractMaster, AbstractMip
export AbstractOracle
export Seq, SeqInOut
export AbstractNorm, StandardNorm, LpNorm
export DisjunctiveCutsAppendRule, NoDisjunctiveCuts, AllDisjunctiveCuts, DisjunctiveCutsSmallerIndices
export SplitIndexSelectionRule, RandomFractional, MostFractional, LargestFractional
export TerminationStatus, NotSolved, TimeLimit, Optimal, InfeasibleOrNumericalIssue
export TimeLimitException, UnexpectedModelStatusException, UndefError
export Hyperplane, aggregate


abstract type AbstractData end
abstract type AbstractBendersDecomposition end
abstract type AbstractMip end
abstract type AbstractMaster end

"""
Subtypes should implement `generate_cuts` to return separating hyperplanes based on the given candidate solutions.
"""
abstract type AbstractOracle end


# ============================================================================
# Global data type; Problem Data is optional; user can define their own structure for problem-specific data
# ============================================================================
struct Data
    dim_x::Int
    dim_t::Int
    problem::AbstractData
    c_x::Vector{Float64}
    c_t::Vector{Float64}
end

# ============================================================================
# Algorithm Strategy Types
# ============================================================================
# # Loop strategies
# struct Sequential <: SolutionProcedure end
# struct Callback <: SolutionProcedure end
# struct StochasticSequential <: SolutionProcedure end
# struct StochasticCallback <: SolutionProcedure end

# # Cut strategies
# struct ClassicalCut <: CutStrategy end
# struct FatKnapsackCut <: CutStrategy end
# struct SlimKnapsackCut <: CutStrategy end
# struct KnapsackCut <: CutStrategy end

# Oracle types
# struct ClassicalOracle <: OracleType end
# struct UserTypicalOracle <: OracleType end
# Loop strategies
abstract type LoopStrategy end
struct Seq <: LoopStrategy end
struct SeqInOut <: LoopStrategy end


# # ============================================================================
# # Cut Strengthening Types
# # ============================================================================
# abstract type CutStrengtheningType end
# struct PureDisjunctiveCut <: CutStrengtheningType end
# struct StrengthenedDisjunctiveCut <: CutStrengtheningType end

# ============================================================================
# Normalization
# ============================================================================
abstract type AbstractNorm end
struct StandardNorm <: AbstractNorm end
mutable struct LpNorm <: AbstractNorm 
    p::Float64
    function LpNorm(p::Float64)
        new(p)
    end
end

# struct L1Norm <: LpNorm end
# struct L2Norm <: LpNorm end
# struct LInfNorm <: LpNorm end

# ============================================================================
# Norm Types
# ============================================================================
abstract type SplitIndexSelectionRule end
abstract type SimpleSplit <: SplitIndexSelectionRule end
struct RandomFractional <: SimpleSplit end
struct MostFractional <: SimpleSplit end
struct LargestFractional <: SimpleSplit end

# # ============================================================================
# # Disjunction System
# # ============================================================================
# struct DisjunctiveCut <: CutStrategy
#     base_cut_strategy::CutStrategy
#     norm_type::AbstractNorm
#     cut_strengthening_type::CutStrengtheningType
#     use_two_sided_cuts::Bool
#     include_master_cuts::Bool
#     reuse_dcglp::Bool
#     verbose::Bool
# end

abstract type TerminationStatus end
struct NotSolved <: TerminationStatus end
struct TimeLimit <: TerminationStatus end
struct Optimal <: TerminationStatus end
struct InfeasibleOrNumericalIssue <: TerminationStatus end

struct TimeLimitException <: Exception 
    msg::String
end

struct UnexpectedModelStatusException <: Exception 
    msg::String
end

struct AlgorithmException <: Exception 
    msg::String
end

struct UndefError <: Exception 
    msg::String
end

abstract type DisjunctiveCutsAppendRule end
struct NoDisjunctiveCuts <: DisjunctiveCutsAppendRule end
struct AllDisjunctiveCuts <: DisjunctiveCutsAppendRule end
struct DisjunctiveCutsSmallerIndices <: DisjunctiveCutsAppendRule end



