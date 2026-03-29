const _CPLEX_OPTIMIZER_FACTORY = Ref{Union{Nothing,Function}}(nothing)
const _GUROBI_OPTIMIZER_FACTORY = Ref{Union{Nothing,Function}}(nothing)

_normalize_optimizer_attribute(key::AbstractString, value) =
    MOI.RawOptimizerAttribute(String(key)) => value
_normalize_optimizer_attribute(key, value) = key => value

function _normalize_optimizer_attributes(args::Pair...)
    return map(arg -> _normalize_optimizer_attribute(arg.first, arg.second), args)
end

function _register_cplex_optimizer!(factory::Function)
    _CPLEX_OPTIMIZER_FACTORY[] = factory
    return nothing
end

function _register_gurobi_optimizer!(factory::Function)
    _GUROBI_OPTIMIZER_FACTORY[] = factory
    return nothing
end

"""
    cplex_optimizer(args::Pair...)

Build a CPLEX optimizer with optional JuMP/MOI attributes.

This helper is provided by the optional `BendersXCPLEXExt` package extension.
Install and load `CPLEX` before calling it.
"""
function cplex_optimizer(args::Pair...)
    factory = _CPLEX_OPTIMIZER_FACTORY[]
    if !isnothing(factory)
        return Base.invokelatest(factory, args...)
    end
    throw(ArgumentError(
        "CPLEX support is optional. Install and load CPLEX first, then call `BendersX.cplex_optimizer(...)`."
    ))
end

"""
    gurobi_optimizer(args::Pair...)

Build a Gurobi optimizer with optional JuMP/MOI attributes.

This helper is provided by the optional `BendersXGurobiExt` package extension.
Install and load `Gurobi` before calling it.
"""
function gurobi_optimizer(args::Pair...)
    factory = _GUROBI_OPTIMIZER_FACTORY[]
    if !isnothing(factory)
        return Base.invokelatest(factory, args...)
    end
    throw(ArgumentError(
        "Gurobi support is optional. Install and load Gurobi first, then call `BendersX.gurobi_optimizer(...)`."
    ))
end

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

"""
    assign_attributes!(model::Model, config::AbstractDict{<:AbstractString,<:Any})

Assign a solver and solver attributes to `model` from a configuration
dictionary.

The dictionary must contain a `"solver"` entry with value `"Gurobi"` or
`"CPLEX"` (case-insensitive). All remaining entries are forwarded as optimizer
attributes. The helper finishes by enabling silent solver output.
"""
function assign_attributes!(model::Model, config::AbstractDict{<:AbstractString,<:Any})
    haskey(config, "solver") || throw(ArgumentError("The configuration dictionary must contain a \"solver\" entry."))

    solver = lowercase(strip(string(config["solver"])))
    optimizer_attributes = Pair{Any,Any}[]
    for (param, value) in config
        if param != "solver"
            push!(optimizer_attributes, param => value)
        end
    end

    if solver == "gurobi"
        set_optimizer(model, gurobi_optimizer(optimizer_attributes...))
    elseif solver == "cplex"
        set_optimizer(model, cplex_optimizer(optimizer_attributes...))
    else
        throw(ArgumentError("Unsupported solver: $(config["solver"])"))
    end

    set_optimizer_attribute(model, MOI.Silent(), true)
    return model
end
