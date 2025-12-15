
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



