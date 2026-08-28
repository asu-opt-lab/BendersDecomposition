# -----------------------------------------------------------------------------
# Lp-distance normalization
# -----------------------------------------------------------------------------

"""
    LpDistanceNormalization <: AbstractNormalization

``L_p``-distance normalization for DCGLP-based disjunctive cut generation.

# Fields

- `norm_p::Float64`: Norm used for normalization. Supported values are `1`, `2`, and `Inf`.

# Constructor

    LpDistanceNormalization(p = Inf)

Construct an ``L_p``-distance normalization using the specified norm.

# Throws

Throws an `ArgumentError` if `p` is not one of `1`, `2`, or `Inf`.

See also: [`AbstractNormalization`](@ref), [`ReversePolarNormalization`](@ref)
"""
mutable struct LpDistanceNormalization <: AbstractNormalization
    norm_p::Float64

    function LpDistanceNormalization(p::Real = Inf)
        p = Float64(p)
        p in (1.0, 2.0, Inf) || throw(ArgumentError("LpDistanceNormalization: Unsupported norm $p"))
        return new(p)
    end
end

function add_normalization_constraint!(
    normalization::LpDistanceNormalization,
    master::AbstractMaster,
    dcglp::Model,
    tau::VariableRef,
    sx::AbstractVector{VariableRef},
    st::AbstractVector{VariableRef},
)

    var_vec = [tau; sx; st]
    p = normalization.norm_p
    if p == 1.0
        @constraint(dcglp, concone, var_vec in MOI.NormOneCone(length(var_vec)))
    elseif p == 2.0
        @constraint(dcglp, concone, var_vec in MOI.SecondOrderCone(length(var_vec)))
    elseif p == Inf
        @constraint(dcglp, concone, var_vec in MOI.NormInfinityCone(length(var_vec)))
    else
        throw(ArgumentError("LpDistanceNormalization: Unsupported norm $p"))
    end
    return nothing
end

function update_dcglp_for_candidate!(normalization::LpDistanceNormalization, dcglp::Model, x_value::Vector{Float64}, t_value::Vector{Float64})
    set_normalized_rhs.(dcglp[:conx], x_value)
    set_normalized_rhs.(dcglp[:cont], t_value)
end

function update_dcglp_upper_bound_and_gap!(
    normalization::LpDistanceNormalization,
    state::DcglpState,
    log::DcglpLog,
    t_value::Vector{Float64},
)
    update_upper_bound_and_gap!(
        state,
        log,
        (t1, t2) -> LinearAlgebra.norm([state.values[:sx]; t1 .+ t2 .- t_value], normalization.norm_p),
    )
    return nothing
end

function disjunctive_cut_normalization_value(
    normalization::LpDistanceNormalization,
    dcglp::Model,
    gamma_x::Vector{Float64},
    gamma_t::Vector{Float64},
)
    p = normalization.norm_p
    if p == 1.0
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), Inf)
    elseif p == 2.0
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), 2.0)
    elseif p == Inf
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), 1.0)
    else
        throw(ArgumentError("LpDistanceNormalization: Unsupported norm $p"))
    end
    return max(1.0, norm_value)
end