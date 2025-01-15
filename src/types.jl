# ============================================================================
# Exports
# ============================================================================
export CFLPData, UFLPData, SCFLPData, SNIPData
export Sequential, Callback, StochasticSequential
export ClassicalCut, FatKnapsackCut, SlimKnapsackCut, KnapsackCut
export PureDisjunctiveCut, StrengthenedDisjunctiveCut
export StandardNorm, L1Norm, L2Norm, LInfNorm
export DisjunctiveCut

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

struct SCFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    n_scenarios::Int
    capacities::Vector{Float64}
    demands::Vector{Vector{Float64}}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

struct MCNDPData <: AbstractData
    num_nodes::Int      # Number of nodes
    num_arcs::Int       # Number of arcs
    num_commodities::Int # Number of commodities
    arcs::Vector{Tuple{Int,Int}}  # Arcs (from_node, to_node)
    fixed_costs::Vector{Float64}   # Fixed costs for opening arcs
    variable_costs::Vector{Float64} # Variable costs per unit flow
    capacities::Vector{Float64}     # Arc capacities
    demands::Vector{Tuple{Int,Int,Float64}} # Demands (origin, destination, quantity)
end

struct SNIPData <: AbstractData
    num_nodes::Int
    num_scenarios::Int
    scenarios::Vector{Tuple{Int,Int,Float64}} # (from_node, to_node, probability)
    D::Vector{Tuple{Int,Int,Float64,Float64}} # (from_node, to_node, r, q)
    A_minus_D::Vector{Tuple{Int,Int,Float64}} # (from_node, to_node, r)
    ψ::Vector{Vector{Float64}} # ψ matrix
    budget::Float64
end

# ============================================================================
# Algorithm Strategy Types
# ============================================================================
# Loop strategies
struct Sequential <: SolutionProcedure end
struct Callback <: SolutionProcedure end
struct StochasticSequential <: SolutionProcedure end
struct StochasticCallback <: SolutionProcedure end

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









