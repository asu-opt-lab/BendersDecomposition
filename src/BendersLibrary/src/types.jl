

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
# Rules for constructing a split set # To-Do: change SplitIndexSelectionRule to SplitSelectionRule
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




