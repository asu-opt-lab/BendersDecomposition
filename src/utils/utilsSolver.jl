# Q. `DEFAULT_OPTIMIZER` and `default_optimizer()` are defining the same thing. why is it defined twice?
const DEFAULT_OPTIMIZER = optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => true)

default_optimizer() = optimizer_with_attributes(GLPK.Optimizer, MOI.Silent() => true)

function set_optimizer_checked!(model::Model, optimizer, context::AbstractString)
    try
        set_optimizer(model, optimizer)
    catch err
        throw(ArgumentError(
            "$(context) could not attach optimizer $(repr(optimizer)). " *
            "Omit `optimizer` to use the default GLPK optimizer, or pass a JuMP-compatible optimizer constructor. " *
            "Original error: $(sprint(showerror, err))"
        ))
    end
    return model
end

callback_solver_tag(model::Model) = Val(Symbol(solver_name(model)))

"""
    callback_node_count(cb_data, model::Model)

Return the solver-reported callback node count, or `nothing` if this metadata is not supported for the active solver.

Solver-specific support is provided by specializing `callback_node_count(cb_data, model, ::Val{:SolverName})` for the solver of interest through a package extension. For instance,CPLEX support is implemented by `BendersXCPLEXExt` when `CPLEX.jl` is loaded. If no implementation is available, this function emits a warning and returns `nothing`.

See also: [`callback_node_depth`](@ref)
"""
function callback_node_count(cb_data, model::Model)
    node_count = callback_node_count(cb_data, model, callback_solver_tag(model))
    if isnothing(node_count)
        message =
            "$(solver_name(model)) callback node count metadata is not currently supported. " *
            "If you need it, add a solver extension implementing `callback_node_count`."
        @warn message maxlog = 1 _id = :callback_missing_node_count_metadata
    end
    return node_count
end
callback_node_count(cb_data, model::Model, ::Val) = nothing

"""
    callback_node_depth(cb_data, model::Model)

Return the solver-reported node depth for the current callback, or nothing if node-depth metadata is not supported for the active solver.

Solver-specific support is provided by specializing `callback_node_depth(cb_data, model, ::Val{:SolverName})` for the solver of interest through a package extension. For instance, CPLEX support is implemented by `BendersXCPLEXExt` when `CPLEX.jl` is loaded. If no implementation is available, this function emits a warning and returns `nothing`.

# Notes
- Solver-specific support is provided through package extensions.
- CPLEX support is implemented by `BendersXCPLEXExt` when `CPLEX.jl` is loaded.
- If the active solver does not support this metadata, BendersX emits a warning
  and returns `nothing`.

See also: [`callback_node_count`](@ref)
"""
function callback_node_depth(cb_data, model::Model)
    node_depth = callback_node_depth(cb_data, model, callback_solver_tag(model))
    if isnothing(node_depth)
        message =
            "$(solver_name(model)) callback node depth metadata is not currently supported. " *
            "If you need it, add a solver extension implementing `callback_node_depth`."
        @warn message maxlog = 1 _id = :callback_missing_depth_metadata
    end
    return node_depth
end
callback_node_depth(cb_data, model::Model, ::Val) = nothing
