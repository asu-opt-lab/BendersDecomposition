function add_normalization_constraint!(
    ::DirectionalReversePolarNormalization,
    dcglp::Model,
    tau::VariableRef,
    sx::AbstractVector{VariableRef},
    st::AbstractVector{VariableRef},
)
    @constraint(dcglp, con_reverse_polar_x[j in eachindex(sx)], sx[j] + 0.0 * tau == 0.0)
    @constraint(dcglp, con_reverse_polar_t[j in eachindex(st)], st[j] + 0.0 * tau == 0.0)
end

function initialize_dcglp_state(::DirectionalReversePolarNormalization)
    state = DcglpState()
    state.values[:oracle_statuses] = fill(:inactive, 2)
    state.values[:oracle_gaps] = fill(0.0, 2)
    state.values[:oracle_scaled_gaps] = fill(0.0, 2)
    return state
end

function should_fallback_typical(normalization::DirectionalReversePolarNormalization, oracle::SplitOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    direction_x = x_value .- normalization.core_point_x
    direction_t = t_value .- normalization.core_point_t

    if LinearAlgebra.norm(vcat(direction_x, direction_t), Inf) <= oracle.param.zero_tol
        return true
    end

    normalization.last_direction_x = direction_x
    normalization.last_direction_t = direction_t
    return false
end

function update_dcglp_for_candidate!(normalization::DirectionalReversePolarNormalization, oracle::SplitOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    set_normalized_rhs.(oracle.dcglp[:cont], t_value)
    update_reverse_polar_constraints!(oracle, normalization.last_direction_x, normalization.last_direction_t)
end

function record_dcglp_oracle_result!(::DirectionalReversePolarNormalization, state::DcglpState, i::Int, t_block::Vector{Float64})
    state.values[:oracle_statuses][i] = state.is_in_L[i] ? :feasible : :infeasible
    state.values[:oracle_gaps][i], state.values[:oracle_scaled_gaps][i] = compute_directional_oracle_gap(
        state.is_in_L[i],
        t_block,
        state.f_x[i],
        state.values[:ω_0][i],
    )
end

function build_dcglp_disjunctive_cut(
    normalization::DirectionalReversePolarNormalization,
    dcglp::Model,
    common::SplitOracleParam,
    ::Float64,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int},
)
    isempty(normalization.last_direction_x) &&
        throw(AlgorithmException("DirectionalReversePolarNormalization: no direction vector cached. Call update_dcglp_for_candidate! first."))

    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    direction_value = normalize_directional_duals!(
        gamma_x,
        gamma_t,
        normalization.last_direction_x,
        normalization.last_direction_t;
        zero_tol = common.zero_tol,
    )

    gamma_0 = dual(dcglp[:con0]) / direction_value
    gamma_x, gamma_0 = apply_lift_or_strengthen(
        dcglp, gamma_x, zero_indices, one_indices;
        lift = common.lift, strengthen = common.strengthened,
        zero_tol = common.zero_tol, gamma_0 = gamma_0,
    )
    return Hyperplane(gamma_x, gamma_t, gamma_0)
end

"""
    set_core_point!(oracle::SplitOracle, core_point_x, core_point_t)

Update the core point used by a directional reverse-polar `SplitOracle`. The supplied
vectors must have the same dimensions as the oracle's existing core point.
"""
function set_core_point!(oracle::SplitOracle, core_point_x::Vector{Float64}, core_point_t::Vector{Float64})
    return set_core_point!(oracle.param.normalization, core_point_x, core_point_t)
end

function set_core_point!(normalization::DirectionalReversePolarNormalization, core_point_x::Vector{Float64}, core_point_t::Vector{Float64})
    length(core_point_x) == length(normalization.core_point_x) ||
        throw(DimensionMismatch("`core_point_x` has length $(length(core_point_x)) but expected $(length(normalization.core_point_x))."))
    length(core_point_t) == length(normalization.core_point_t) ||
        throw(DimensionMismatch("`core_point_t` has length $(length(core_point_t)) but expected $(length(normalization.core_point_t))."))

    normalization.core_point_x .= core_point_x
    normalization.core_point_t .= core_point_t
end

function update_reverse_polar_constraints!(
    oracle::SplitOracle,
    direction_x::Vector{Float64},
    direction_t::Vector{Float64},
)
    dcglp = oracle.dcglp

    if haskey(dcglp, :con_reverse_polar_x) && haskey(dcglp, :con_reverse_polar_t) &&
       length(dcglp[:con_reverse_polar_x]) == length(direction_x) && length(dcglp[:con_reverse_polar_t]) == length(direction_t)
        tau = dcglp[:tau]
        for j in eachindex(direction_x)
            set_normalized_coefficient(dcglp[:con_reverse_polar_x][j], tau, direction_x[j])
        end
        for j in eachindex(direction_t)
            set_normalized_coefficient(dcglp[:con_reverse_polar_t][j], tau, direction_t[j])
        end
        return nothing
    end

    delete_registered_constraints!(dcglp, :con_reverse_polar_x)
    delete_registered_constraints!(dcglp, :con_reverse_polar_t)

    @constraint(
        dcglp,
        con_reverse_polar_x[j in eachindex(direction_x)],
        dcglp[:sx][j] + direction_x[j] * dcglp[:tau] == 0.0,
    )
    @constraint(
        dcglp,
        con_reverse_polar_t[j in eachindex(direction_t)],
        dcglp[:st][j] + direction_t[j] * dcglp[:tau] == 0.0,
    )
end

function compute_directional_oracle_gap(
    is_in_L::Bool,
    t_block::Vector{Float64},
    f_x_i::Vector{Float64},
    omega_0::Float64,
)
    if is_in_L
        return 0.0, 0.0
    elseif any(isnan, f_x_i)
        return NaN, NaN
    end

    gap = maximum(max.(f_x_i .- t_block, 0.0))
    return gap, omega_0 * gap
end
