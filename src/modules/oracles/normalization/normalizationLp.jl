# -----------------------------------------------------------------------------
# Lp-distance normalization
# -----------------------------------------------------------------------------

"""
    LpDistanceNormalization(p)

Distance-norm DCGLP normalization for `SplitOracle`.
"""
mutable struct LpDistanceNormalization <: AbstractDisjunctiveNormalization
    norm_p::Float64

    function LpDistanceNormalization(p::Real = Inf)
        p = Float64(p)
        p in (1.0, 2.0, Inf) || throw(ArgumentError("Unsupported norm p: $p"))
        return new(p)
    end
end

function add_normalization_constraint!(
    normalization::LpDistanceNormalization,
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
        throw(UndefError("Unsupported norm p: $p"))
    end
end

function update_dcglp_upper_bound_and_gap!(
    normalization::LpDistanceNormalization,
    state::DcglpState,
    log::DcglpLog,
    reference_t::Vector{Float64},
    t_value::Vector{Float64},
)
    fill_dcglp_omega_t_estimates!(state, t_value)
    all(f_i -> !any(isnan, f_i), state.f_x) || return nothing
    update_upper_bound_and_gap!(
        state,
        log,
        (t1, t2) -> LinearAlgebra.norm([state.values[:sx]; t1 .+ t2 .- reference_t], normalization.norm_p),
    )
end

function lp_cut_dual_norm_value(gamma_x, gamma_t, p::Float64)
    if p == 1.0
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), Inf)
    elseif p == 2.0
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), 2.0)
    elseif p == Inf
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), 1.0)
    else
        throw(UndefError("Unsupported norm p: $p"))
    end
    return max(1.0, norm_value)
end

function disjunctive_cut_normalization_value(
    normalization::LpDistanceNormalization,
    gamma_x::Vector{Float64},
    gamma_t::Vector{Float64},
    common::SplitOracleParam,
    current_lb::Float64,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int},
)
    common.lift && (!isempty(zero_indices) || !isempty(one_indices)) ||
        return 1.0
    return lp_cut_dual_norm_value(gamma_x, gamma_t, normalization.norm_p)
end
