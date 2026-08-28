# ============================================================================
# Benders execution environments
# ============================================================================

# ---------------------------------------------------------------------------- 
# Common interface 
# ----------------------------------------------------------------------------
"""
    solve!(env::AbstractBendersEnv; kwargs...)

Execute a Benders decomposition environment.

Concrete execution environments implement this method according to their algorithmic workflow.

# Error Handling

The execution environment handles the following internal exceptions:

- [`TimeLimitException`](@ref): Sets `env.termination_status` to [`TimeLimit`](@ref).
- [`UnexpectedModelStatusException`](@ref): Sets `env.termination_status` to [`InfeasibleOrNumericalIssue`](@ref).

These exceptions may originate directly from `solve!` or propagate from operations invoked during the solve, such as preprocessing or oracle evaluation.

Other exceptions are not converted to algorithm-level termination statuses and propagate to the caller.

# Returns

The return value is defined by the concrete execution environment.

See also: [`AbstractBendersEnv`](@ref), [`TerminationStatus`](@ref)
"""
function solve!(env::AbstractBendersEnv)
    throw(UnimplementedInterfaceException("update solve! for $(typeof(env))"))
end

# ---------------------------------------------------------------------------- 
# Benders LP relaxation preprocessing
# ----------------------------------------------------------------------------
include("preprocessing/preprocessing.jl") # must be included first

# ---------------------------------------------------------------------------- 
# Concrete environment implementations 
# ----------------------------------------------------------------------------
include("BendersSeq.jl")
include("BendersBnB.jl")
include("BendersSeqInOut.jl")