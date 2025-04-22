export EmptyCallbackParam, lazy_callback, UserCallback, NoUserCallback

"""
    AbstractCallbackParam

Abstract type for parameters used in callbacks during the branch-and-bound process.
These parameters control how and when callbacks are executed.
"""
abstract type AbstractCallbackParam end

"""
    EmptyCallbackParam <: AbstractCallbackParam

Represents empty (default) parameters for callbacks.
Used when no specific parameters are needed for a callback.
"""
struct EmptyCallbackParam <: AbstractCallbackParam
end

"""
    AbstractLazyCallback

Abstract type for lazy constraint callbacks in Benders decomposition.
Lazy callbacks are used to dynamically add cuts when integer solutions are found
during the branch-and-bound process.
"""
abstract type AbstractLazyCallback end

"""
    lazy_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::AbstractLazyCallback)

Generic function for implementing lazy callbacks. This should be overridden by specific implementations of `AbstractLazyCallback`.

# Arguments
- `cb_data`: Callback data from the solver
- `master_model::Model`: The JuMP master problem model
- `log::BendersBnBLog`: Log object to record statistics
- `callback::AbstractLazyCallback`: Configuration for the lazy callback
"""
function lazy_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::AbstractLazyCallback)
    throw(UndefError("Lazy callback not implemented for $(typeof(callback))"))
end


"""
    AbstractUserCallback

Abstract type for user cut callbacks in Benders decomposition.
User callbacks are used to dynamically add cuts at fractional nodes during the branch-and-bound process to strengthen the formulation.
"""
abstract type AbstractUserCallback end

"""
    NoUserCallback <: AbstractUserCallback

Represents a no-operation user callback.
Use this when you don't want to add any user cuts during the branch-and-bound process.
"""
struct NoUserCallback <: AbstractUserCallback end

"""
    user_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::AbstractUserCallback)

Generic function for implementing user callbacks. This should be overridden by specific implementations of `AbstractUserCallback`.

# Arguments
- `cb_data`: Callback data from the solver
- `master_model::Model`: The JuMP master problem model
- `log::BendersBnBLog`: Log object to record statistics
- `callback::AbstractUserCallback`: Configuration for the user callback
"""
function user_callback(cb_data, master_model::Model, log::BendersBnBLog, callback::AbstractUserCallback)
    # Default implementation for the abstract type
    # Concrete subtypes should override this method except for NoUserCallback
    if callback isa NoUserCallback
        return  # Silent no-op for NoUserCallback
    end
    throw(UndefError("User callback not implemented for $(typeof(callback))"))
end

include("callbackLazy.jl")
include("callbackUser.jl")

