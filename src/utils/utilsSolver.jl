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

function callback_node_count(cb_data, model::Model)
    return callback_node_count(cb_data, model, callback_solver_tag(model))
end
callback_node_count(cb_data, model::Model, ::Val) = nothing

function callback_node_depth(cb_data, model::Model)
    return callback_node_depth(cb_data, model, callback_solver_tag(model))
end
callback_node_depth(cb_data, model::Model, ::Val) = nothing
