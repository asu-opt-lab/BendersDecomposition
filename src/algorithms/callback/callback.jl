export EmptyCallbackParam

"""
    AbstractCallbackParam

Abstract type for parameters used in callbacks during the branch-and-bound process.
"""
abstract type AbstractCallbackParam end

"""
    EmptyCallbackParam <: AbstractCallbackParam

Represents empty (default) parameters for callbacks.
"""
struct EmptyCallbackParam <: AbstractCallbackParam
end

"""
    AbstractLazyCallback

Abstract type for lazy constraint callbacks in Benders decomposition.
"""
abstract type AbstractLazyCallback end

"""
    AbstractUserCallback

Abstract type for user cut callbacks in Benders decomposition.
"""
abstract type AbstractUserCallback end

include("callbackLazy.jl")
include("callbackUser.jl")

