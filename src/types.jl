# ============================================================================
# Abstract type hierarchy
# ============================================================================
abstract type AbstractData end
abstract type SolutionProcedure end
abstract type CutStrategy end

abstract type AbstractMasterProblem end
abstract type AbstractSubProblem end
abstract type AbstractDCGLP end
abstract type AbstractMILP end

# ============================================================================
# Problem Data Types
# ============================================================================
struct CFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    capacities::Vector{Float64}
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

struct UFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

# ============================================================================
# Algorithm Strategy Types
# ============================================================================
# Loop strategies
struct Sequential <: SolutionProcedure end
struct Callback <: SolutionProcedure end

# Cut strategies
struct ClassicalCut <: CutStrategy end
struct FatKnapsackCut <: CutStrategy end
struct SlimKnapsackCut <: CutStrategy end
struct KnapsackCut <: CutStrategy end

# ============================================================================
# Cut Strengthening Types
# ============================================================================
abstract type CutStrengtheningType end
struct PureDisjunctiveCut <: CutStrengtheningType end
struct StrengthenedDisjunctiveCut <: CutStrengtheningType end

# ============================================================================
# Norm Types
# ============================================================================
abstract type AbstractNormType end
struct StandardNorm <: AbstractNormType end
abstract type LNorm <: AbstractNormType end

struct L1Norm <: LNorm end
struct L2Norm <: LNorm end
struct LInfNorm <: LNorm end

# ============================================================================
# Disjunction System
# ============================================================================
struct DisjunctiveCut <: CutStrategy
    base_cut_strategy::CutStrategy
    norm_type::AbstractNormType
    cut_strengthening_type::CutStrengtheningType
    use_two_sided_cuts::Bool
    include_master_cuts::Bool
    reuse_dcglp::Bool
    verbose::Bool
end

# ============================================================================
# Exports
# ============================================================================
export CFLPData, UFLPData
export Sequential, Callback
export ClassicalCut, FatKnapsackCut, SlimKnapsackCut, KnapsackCut
export PureDisjunctiveCut, StrengthenedDisjunctiveCut
export StandardNorm, L1Norm, L2Norm, LInfNorm
export DisjunctiveCut







