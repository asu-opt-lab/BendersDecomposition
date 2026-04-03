function attach_optimizer!(model::Model, optimizer)
    if !isnothing(optimizer)
        set_optimizer(model, optimizer)
    end
    return model
end

function has_optimizer_attached(model::Model)
    try
        MOI.get(backend(model), MOI.SolverName())
        return true
    catch err
        if err isa MOI.GetAttributeNotAllowed
            return false
        end
        rethrow()
    end
end

function require_optimizer_attached(model::Model, context::AbstractString)
    if !has_optimizer_attached(model)
        throw(ArgumentError(
            "$(context) requires an optimizer, but none is attached. " *
            "Pass `optimizer = ...` when constructing `Master` or model-based oracles, " *
            "or attach an optimizer to the `Model` before calling `customize_mip_model!` / `optimize!`."
        ))
    end
    return model
end

callback_solver_tag(model::Model) = Val(Symbol(solver_name(model)))

# Solver extensions override these methods to expose callback metadata.
# Returning `nothing` means the active solver does not provide the requested
# metadata in the callback context.
function callback_node_count(cb_data, model::Model)
    node_count = callback_node_count(cb_data, model, callback_solver_tag(model))
    if isnothing(node_count)
        message =
            "$(solver_name(model)) does not provide callback node count metadata. " *
            "If this solver supports node count metadata, add a solver extension implementing `callback_node_count`."
        @warn message maxlog = 1 _id = :callback_missing_node_count_metadata
    end
    return node_count
end
callback_node_count(cb_data, model::Model, ::Val) = nothing

function callback_node_depth(cb_data, model::Model)
    node_depth = callback_node_depth(cb_data, model, callback_solver_tag(model))
    if isnothing(node_depth)
        message =
            "$(solver_name(model)) does not provide callback node depth metadata. " *
            "If this solver supports node depth metadata, add a solver extension implementing `callback_node_depth`."
        @warn message maxlog = 1 _id = :callback_missing_depth_metadata
    end
    return node_depth
end
callback_node_depth(cb_data, model::Model, ::Val) = nothing
