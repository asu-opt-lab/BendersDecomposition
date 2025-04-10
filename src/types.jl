# ============================================================================
# Exports
# ============================================================================
# export CFLPData, UFLPData, SCFLPData, MCNDPData,SNIPData
# export Sequential, Callback, StochasticSequential, StochasticCallback
# export ClassicalCut, FatKnapsackCut, SlimKnapsackCut, KnapsackCut
# export PureDisjunctiveCut, StrengthenedDisjunctiveCut
export StandardNorm, LpNorm, AbstractNorm
# export DisjunctiveCut

export Data, AbstractData
export AbstractMaster, AbstractMip
export AbstractOracle, AbstractTypicalOracle, AbstractDisjunctiveOracle
export Seq, SeqInOut
export RandomFractional, MostFractional, LargestFractional
export TerminationStatus, NotSolved, TimeLimit, Optimal, InfeasibleOrNumericalIssue
export TimeLimitException, UnexpectedModelStatusException
# ============================================================================
# Abstract type hierarchy
# ============================================================================
# abstract type SolutionProcedure end
# abstract type CutStrategy end

# abstract type AbstractMasterProblem end
# abstract type AbstractSubProblem end
# abstract type AbstractDCGLP end
# abstract type AbstractMILP end

abstract type AbstractData end
abstract type OracleType end
abstract type LoopStrategy end
abstract type AbstractMip end
abstract type AbstractMaster end
abstract type AbstractOracle end
abstract type AbstractTypicalOracle <: AbstractOracle end
abstract type AbstractDisjunctiveOracle <: AbstractOracle end


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
# Problem Data Types
# ============================================================================
# struct CFLPData <: AbstractData
#     n_facilities::Int
#     n_customers::Int
#     capacities::Vector{Float64}
#     demands::Vector{Float64}
#     fixed_costs::Vector{Float64}
#     costs::Matrix{Float64}
# end

# struct UFLPData <: AbstractData
#     n_facilities::Int
#     n_customers::Int
#     demands::Vector{Float64}
#     fixed_costs::Vector{Float64}
#     costs::Matrix{Float64}
# end

# struct SCFLPData <: AbstractData
#     n_facilities::Int
#     n_customers::Int
#     num_scenarios::Int
#     capacities::Vector{Float64}
#     demands::Vector{Vector{Float64}}
#     fixed_costs::Vector{Float64}
#     costs::Matrix{Float64}
# end

# struct MCNDPData <: AbstractData
#     num_nodes::Int      # Number of nodes
#     num_arcs::Int       # Number of arcs
#     num_commodities::Int # Number of commodities
#     arcs::Vector{Tuple{Int,Int}}  # Arcs (from_node, to_node)
#     fixed_costs::Vector{Float64}   # Fixed costs for opening arcs
#     variable_costs::Vector{Float64} # Variable costs per unit flow
#     capacities::Vector{Float64}     # Arc capacities
#     demands::Vector{Tuple{Int,Int,Float64}} # Demands (origin, destination, quantity)
# end

# struct SNIPData <: AbstractData
#     num_nodes::Int
#     num_scenarios::Int
#     scenarios::Vector{Tuple{Int,Int,Float64}} # (from_node, to_node, probability)
#     D::Vector{Tuple{Int,Int,Float64,Float64}} # (from_node, to_node, r, q)
#     A_minus_D::Vector{Tuple{Int,Int,Float64}} # (from_node, to_node, r)
#     ψ::Vector{Vector{Float64}} # ψ matrix
#     budget::Float64
# end

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
struct RandomFractional <: SplitIndexSelectionRule end
struct MostFractional <: SplitIndexSelectionRule end
struct LargestFractional <: SplitIndexSelectionRule end

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


mutable struct BendersState
    iteration::Int
    master_time::Float64
    oracle_time::Float64
    total_time::Float64
    is_in_L::Bool
    LB::Float64
    UB::Float64
    gap::Float64
   
    # Constructor with specified values
    function BendersState()
        new(0, 0.0, 0.0, 0.0, false, -Inf, Inf, 100.0)
    end
end

mutable struct BendersDecompositionLog
    iterations::Vector{BendersState}
    start_time::Float64
    master_time::Float64
    oracle_time::Float64
    LB::Float64
    UB::Float64
    termination_status::TerminationStatus
    
    function BendersDecompositionLog()
        new(Vector{BendersState}(), time(), 0.0, 0.0, -Inf, Inf, NotSolved())
    end
end




