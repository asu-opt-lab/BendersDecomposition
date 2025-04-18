# ============================================================================
# Exports
# ============================================================================

export AbstractBendersDecomposition
export Data, AbstractData
export AbstractMaster, AbstractMip
export AbstractOracle, AbstractOracleParam
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
Any concrete subtype of `AbstractOracle` must have a field `oracle_param<:AbstractOracleParam` containing adjustable parameters that affect the oracle's behavior. 
The type of `oracle_param` can be EmptyOracleParam when there is no adjustable oracle parameter. 

Subtypes should implement `generate_cuts` to return separating hyperplanes based on the given candidate solutions.
"""
abstract type AbstractOracle end
abstract type AbstractOracleParam end


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



