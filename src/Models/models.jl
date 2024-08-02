
# Abstract types
abstract type AbstractMasterEnv end
abstract type AbstractDCGLPEnv end
abstract type AbstractSubEnv end
abstract type AbstractMipEnv end

mutable struct SplitInfo 
    # indices
    γ₀s
    γₓs
    γₜs
    # ifaddall::Bool
end

include("CFLP/cflp.jl")
include("UFLP/uflp.jl")
