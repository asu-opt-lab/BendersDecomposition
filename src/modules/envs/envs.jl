# ============================================================================
# Benders execution environments
# ============================================================================

# ---------------------------------------------------------------------------- 
# Common interface 
# ----------------------------------------------------------------------------
"""
    solve!(env::AbstractBendersEnv)

Execute the Benders algorithm defined by `env` and return its execution log.

Concrete [`AbstractBendersEnv`](@ref) types implement this method according to their corresponding Benders algorithm. This fallback method is called only when no implementation is available for the concrete environment type.
"""
function solve!(env::AbstractBendersEnv)
    throw(UndefError("update solve! for $(typeof(env))"))
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
include("SpecializedBendersSeq.jl")