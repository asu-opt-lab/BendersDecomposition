function delete_registered_constraints!(model::Model, sym::Symbol)
    haskey(model, sym) || return nothing
    registered = model[sym]
    if registered isa AbstractArray
        delete.(Ref(model), registered)
    else
        delete(model, registered)
    end
    unregister(model, sym)
end


function fallback_typical_or_throw(
    oracle::AbstractDisjunctiveOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
    msg::String;
    throw_typical_cuts_for_errors::Bool,
)
    if throw_typical_cuts_for_errors
        @warn msg
        return generate_cuts(
            oracle.typical_oracles[1],
            x_value,
            t_value;
            time_limit = get_sec_remaining(start_time, time_limit),
        )
    end
    throw(UnexpectedModelStatusException(msg))
end
