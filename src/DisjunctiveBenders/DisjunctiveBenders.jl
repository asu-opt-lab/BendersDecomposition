export DCGLP
export DisjunctiveCut
export ClassicalCut
export L1Norm
export L2Norm
export LInfNorm
export PureDisjunctiveCut
export StrengthenedDisjunctiveCut
export StandardNorm

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
struct DisjunctiveCut <: AbstractCutStrategy
    base_cut_strategy::AbstractCutStrategy
    norm_type::AbstractNormType
    cut_strengthening_type::CutStrengtheningType
    use_two_sided_cuts::Bool
    include_master_cuts::Bool
    reuse_dcglp::Bool
    verbose::Bool
end


mutable struct DCGLP 
    model::Model
    γ_constraints::Dict{Symbol,Any}
    γ_values::Vector{Tuple{Float64, Vector{Float64}, Union{Float64, Vector{Float64}}}}
    disjunctive_inequalities_constraints::Vector{ConstraintRef}
    dcglp_constraints::Vector{ConstraintRef}
    master_cuts::Union{Vector{AffExpr}, AffExpr}
end


include("base_dcglp.jl")
include("modeling.jl")
include("generate_cuts.jl")
include("update_dcglp_utilities.jl")
include("generate_cuts_utilities.jl")
include("merge_cuts_utilities.jl")