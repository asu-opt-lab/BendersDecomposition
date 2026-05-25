const DCGLP_OMEGA_T_LOWER_BOUND = -1.0e6

function normalize_add_benders_cuts_to_master(add_benders_cuts_to_master::Union{Bool, Int})
    return add_benders_cuts_to_master === true ? 1 :
           add_benders_cuts_to_master === false ? 0 :
           add_benders_cuts_to_master in (0, 1, 2) ? add_benders_cuts_to_master :
           throw(ArgumentError("`add_benders_cuts_to_master` must be true, false, or an integer in {0, 1, 2}."))
end

function validate_fraction_of_benders_cuts_to_master(fraction::Float64)
    0.0 < fraction <= 1.0 ||
        throw(ArgumentError("`fraction_of_benders_cuts_to_master` must lie in (0, 1]."))
    return fraction
end

function validate_two_typical_oracles!(typical_oracles, oracle_name::String)
    length(typical_oracles) == 2 ||
        throw(ArgumentError("$oracle_name requires exactly two typical oracles."))
end

function validate_binary_master!(master::AbstractMaster, oracle_name::String)
    for xi in master.x
        is_binary(xi) ||
            throw(ArgumentError("$oracle_name requires all master variables to be binary."))
    end
end

function initialize_disjunctive_cut_storage(master::AbstractMaster)
    return [Vector{Hyperplane}() for _ in 1:master.dim_x],
           Hyperplane[],
           Vector{Tuple{SparseVector{Float64, Int}, Float64}}()
end

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

"""
    build_dcglp_skeleton!(dcglp::Model, master::AbstractMaster)

Create the shared DCGLP variables and constraints that every disjunctive oracle
relies on: `omega_0`, `omega_x`, `omega_t`, the `con0` simplex constraint, the
`coneta`/`condelta` cut, the `omega_t` lower-bound row, and the two
`transfer_scaled_linear_rows_and_bounds_with_types!` blocks.

Caller is responsible for the oracle-specific `tau`, objective, `conx`/`cont`
constraints (which involve RHS), and any normalization constraint.
"""
function build_dcglp_skeleton!(dcglp::Model, master::AbstractMaster)
    @variable(dcglp, omega_0[1:2] >= 0)

    @variable(dcglp, omega_x[1:2, 1:master.dim_x])
    @variable(dcglp, omega_t[1:2, 1:master.dim_t])

    @constraint(dcglp, [i in 1:2], omega_t[i, :] .>= DCGLP_OMEGA_T_LOWER_BOUND .* omega_0[i])
    @constraint(dcglp, coneta[i in 1:2, j in 1:master.dim_x], 0 >= -omega_0[i] + omega_x[i, j])
    @constraint(dcglp, condelta[i in 1:2, j in 1:master.dim_x], 0 >= -omega_x[i, j])

    @constraint(dcglp, con0, omega_0[1] + omega_0[2] == 1)

    for i in 1:2
        transfer_scaled_linear_rows_and_bounds_with_types!(
            master.model,
            master.x,
            dcglp,
            omega_x[i, :],
            omega_0[i],
        )
    end

    return dcglp
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
