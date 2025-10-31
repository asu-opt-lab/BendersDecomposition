
abstract type AbstractBendersDecomposition end
abstract type AbstractMip end
abstract type AbstractMaster end
"""
Any concrete subtype of `AbstractOracle` must have a field `oracle_param<:AbstractOracleParam` containing adjustable parameters that affect the oracle's behavior.
The type of `oracle_param` can be BasicOracleParam when there is no oracle-specific adjustable parameter.

Subtypes should implement `generate_cuts` to return separating hyperplanes based on the given candidate solutions.
"""
abstract type AbstractOracle end
abstract type AbstractOracleParam end


# ============================================================================
# Global data type; Problem Data is optional; user can define their own structure for problem-specific data
# To-Do: think about the type for `problem`. Should we remove `AbstractData`?
# ============================================================================
abstract type AbstractData end
struct Data{T<:AbstractData}
    dim_x::Int
    dim_t::Int
    problem::T
    c_x::Vector{Float64}
    c_t::Vector{Float64}
end

abstract type AbstractBendersSeq <: AbstractBendersDecomposition end
abstract type AbstractBendersCallback <: AbstractBendersDecomposition end


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



