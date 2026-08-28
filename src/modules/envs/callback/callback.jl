
"""
    AbstractLazyCallback

Abstract supertype for lazy-constraint callbacks used by Benders branch-and-bound environments.

A lazy callback is invoked when the MIP solver encounters a candidate integer solution and may add Benders cuts as lazy constraints.
"""
abstract type AbstractLazyCallback end

"""
    lazy_callback(
        cb_data, 
        master_model::Model, 
        log::AbstractBendersBnBLog, 
        param::AbstractBendersBnBParam, 
        callback::AbstractLazyCallback
    )

Apply the lazy-constraint callback defined by `callback`.

Concrete [`AbstractLazyCallback`](@ref) subtypes should implement this method to evaluate the candidate solution and add the corresponding Benders cuts to `master_model`. See, for example,[`LazyCallback`](@ref).

# Arguments
- `cb_data`: Callback data supplied by the MIP solver.
- `master_model::Model`: JuMP master model
- `log::AbstractBendersBnBLog`: Branch-and-bound execution log.
- `param::AbstractBendersBnBParam`: Parameters controlling the Benders branch-and-bound procedure.
- `callback::AbstractLazyCallback`: Lazy-callback configuration.
"""
function lazy_callback(cb_data, master_model::Model, log::AbstractBendersBnBLog, param::AbstractBendersBnBParam, callback::AbstractLazyCallback)
    throw(UnimplementedInterfaceException("lazy_callback is not implemented for $(typeof(callback))"))
end

"""
    AbstractUserCallbackParam

Abstract supertype for parameters controlling user-cut callbacks.
"""
abstract type AbstractUserCallbackParam end

"""
    AbstractUserCallback

Abstract supertype for user-cut callbacks used by Benders branch-and-bound environments.

A user callback may be invoked at fractional branch-and-bound nodes to generate Benders cuts that strengthen the master formulation.
"""
abstract type AbstractUserCallback end

"""
    NoUserCallback <: AbstractUserCallback

No-op user-cut callback.

Use this type when no user cuts are to be generated during the branch-and-bound process.
"""
struct NoUserCallback <: AbstractUserCallback end

"""
    user_callback(
        cb_data, 
        master_model::Model, 
        log::AbstractBendersBnBLog, 
        param::AbstractBendersBnBParam, 
        callback::AbstractUserCallback
    )

Apply the user-cut callback defined by `callback`.

Concrete [`AbstractUserCallback`](@ref) subtypes should implement this method to evaluate the current fractional solution and add the corresponding Benders cuts to `master_model`.

[`NoUserCallback`](@ref) provides a no-op implementation.

# Arguments
- `cb_data`: Callback data supplied by the MIP solver.
- `master_model::Model`: JuMP master model.
- `log::AbstractBendersBnBLog`: Branch-and-bound execution log.
- `param::AbstractBendersBnBParam`: Parameters controlling the Benders branch-and-bound procedure.
- `callback::AbstractUserCallback`: User-callback configuration.
"""
function user_callback(cb_data, master_model::Model, log::AbstractBendersBnBLog, param::AbstractBendersBnBParam, callback::AbstractUserCallback)
    callback isa NoUserCallback && return # Silent no-op for NoUserCallback

    throw(UnimplementedInterfaceException("user_callback is not implemented for $(typeof(callback))"))
end

# ----------------------------------------------------------------------------
# Callback implementations
# ----------------------------------------------------------------------------

include("callbackLazy.jl")
include("callbackUser.jl")

