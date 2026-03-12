# ============================================================================
# GBC (Generalized Bound Constraints) Bound Type
# ============================================================================
@enum GBCBoundType begin
    UpperBound   # y ≤ f(x)
    LowerBound   # y ≥ f(x)
    Fixed   # y = f(x)
end

abstract type AbstractBendersEnv end
abstract type AbstractMaster end
"""
Any concrete subtype of `AbstractOracle` must have a field `param<:AbstractOracleParam` containing adjustable parameters that affect the oracle's behavior.
The type of `param` can be BasicOracleParam when there is no oracle-specific adjustable parameter.

Subtypes should implement `generate_cuts` to return separating hyperplanes based on the given candidate solutions.
"""
abstract type AbstractOracle end
abstract type AbstractOracleParam end
abstract type AbstractData end

abstract type AbstractBendersSeq <: AbstractBendersEnv end
abstract type AbstractBendersBnB <: AbstractBendersEnv end

# ============================================================================
# Termination Status of Benders Decomposition
# ============================================================================
abstract type TerminationStatus end
struct NotSolved <: TerminationStatus end
struct TimeLimit <: TerminationStatus end
struct Optimal <: TerminationStatus end
struct InfeasibleOrNumericalIssue <: TerminationStatus end

# ============================================================================
# Error Exceptions
# ============================================================================
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
struct UnsupportedModelException <: Exception 
    msg::String
end

# ============================================================================
# Normalization type for CGLP
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
# Rules for constructing a split set
# ============================================================================
abstract type SplitIndexSelectionRule end
abstract type SimpleSplit <: SplitIndexSelectionRule end
struct RandomFractional <: SimpleSplit end
struct MostFractional <: SimpleSplit end
struct LargestFractional <: SimpleSplit end

# ============================================================================
# Rules for appending pre-found disjunctive cuts to dcglp
# ============================================================================
abstract type DisjunctiveCutsAppendRule end
struct NoDisjunctiveCuts <: DisjunctiveCutsAppendRule end
struct AllDisjunctiveCuts <: DisjunctiveCutsAppendRule end
struct DisjunctiveCutsSmallerIndices <: DisjunctiveCutsAppendRule end


