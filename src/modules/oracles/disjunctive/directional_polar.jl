function validate_strategy_specific!(strategy::DirectionalPolarStrategy, master::AbstractMaster)
    length(strategy.core_point_x) == master.dim_x ||
        throw(DimensionMismatch("`core_point_x` has length $(length(strategy.core_point_x)) but expected $(master.dim_x)."))
    length(strategy.core_point_t) == master.dim_t ||
        throw(DimensionMismatch("`core_point_t` has length $(length(strategy.core_point_t)) but expected $(master.dim_t)."))
end

function build_strategy_dcglp(strategy::DirectionalPolarStrategy, master::AbstractMaster, param::DcglpOracleParam{DirectionalPolarStrategy})
    dcglp = Model(param.dcglp_param.optimizer)

    @variable(dcglp, tau >= 0.0)

    @objective(dcglp, Min, tau)

    build_dcglp_skeleton!(dcglp, master)

    # The `tau` coefficient on `conx`/`cont` is updated each call via
    # `set_normalized_coefficient!`. Declaring it explicitly here registers
    # the variable in the sparse pattern so coefficient updates are O(1).
    @constraint(
        dcglp,
        conx[j in 1:master.dim_x],
        dcglp[:omega_x][1, j] + dcglp[:omega_x][2, j] + 0.0 * tau == strategy.core_point_x[j],
    )
    @constraint(
        dcglp,
        cont[j in 1:master.dim_t],
        dcglp[:omega_t][1, j] + dcglp[:omega_t][2, j] + 0.0 * tau == strategy.core_point_t[j],
    )

    return dcglp
end

function initialize_dcglp_state(::DirectionalPolarStrategy)
    state = DcglpState()
    state.values[:oracle_statuses] = fill(:inactive, 2)
    state.values[:oracle_gaps] = fill(0.0, 2)
    state.values[:oracle_scaled_gaps] = fill(0.0, 2)
    return state
end

function should_fallback_typical(strategy::DirectionalPolarStrategy, oracle::DirectionalPolarOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    direction_x = x_value .- strategy.core_point_x
    direction_t = t_value .- strategy.core_point_t

    if LinearAlgebra.norm(vcat(direction_x, direction_t), Inf) <= oracle.param.zero_tol
        return true
    end

    strategy.last_direction_x = direction_x
    strategy.last_direction_t = direction_t
    return false
end

function update_dcglp_for_candidate!(strategy::DirectionalPolarStrategy, oracle::DirectionalPolarOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    update_direction_constraints!(oracle, x_value, t_value, strategy.last_direction_x, strategy.last_direction_t)
end

dcglp_tau_value(::DirectionalPolarStrategy, dcglp::Model) = value(dcglp[:tau])

function record_dcglp_oracle_result!(::DirectionalPolarStrategy, state::DcglpState, i::Int, t_block::Vector{Float64})
    state.values[:oracle_statuses][i] = state.is_in_L[i] ? :feasible : :infeasible
    state.values[:oracle_gaps][i], state.values[:oracle_scaled_gaps][i] = compute_directional_oracle_gap(
        state.is_in_L[i],
        t_block,
        state.f_x[i],
        state.values[:ω_0][i],
    )
end

function update_dcglp_upper_bound_and_gap!(
    ::DirectionalPolarStrategy,
    state::DcglpState,
    log::DcglpLog,
    ::Vector{Float64},
    t_value::Vector{Float64},
)
    fill_dcglp_omega_t_estimates!(state, t_value)
    previous_ub = log.n_iter >= 1 ? log.iterations[end].UB : 1.0
    state.UB = state.is_in_L[1] && state.is_in_L[2] ? min(previous_ub, state.LB) : previous_ub
    state.gap =
        isfinite(state.UB) ?
        max(state.UB - state.LB, 0.0) / max(abs(state.UB), 1.0) * 100 :
        Inf
end

function has_dcglp_disjunctive_cut(::DirectionalPolarStrategy, current_lb::Float64, ::Vector{Float64}, zero_tol::Float64)
    return current_lb >= zero_tol
end

function build_dcglp_disjunctive_cut(
    strategy::DirectionalPolarStrategy,
    dcglp::Model,
    common::DcglpOracleParam{DirectionalPolarStrategy},
    ::Float64,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int},
)
    isempty(strategy.last_direction_x) &&
        throw(AlgorithmException("DirectionalPolarStrategy: no direction vector cached. Call update_dcglp_for_candidate! first."))

    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    direction_value = normalize_directional_duals!(
        gamma_x,
        gamma_t,
        strategy.last_direction_x,
        strategy.last_direction_t;
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
    set_core_point!(oracle::DirectionalPolarOracle, core_point_x, core_point_t)

Update the core point used by a [`DirectionalPolarOracle`](@ref). The supplied
vectors must have the same dimensions as the oracle's existing core point.
"""
function set_core_point!(oracle::DirectionalPolarOracle, core_point_x::Vector{Float64}, core_point_t::Vector{Float64})
    strategy = oracle.param.strategy
    length(core_point_x) == length(strategy.core_point_x) ||
        throw(DimensionMismatch("`core_point_x` has length $(length(core_point_x)) but expected $(length(strategy.core_point_x))."))
    length(core_point_t) == length(strategy.core_point_t) ||
        throw(DimensionMismatch("`core_point_t` has length $(length(core_point_t)) but expected $(length(strategy.core_point_t))."))

    strategy.core_point_x .= core_point_x
    strategy.core_point_t .= core_point_t
end

function update_direction_constraints!(
    oracle::DirectionalPolarOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    direction_x::Vector{Float64},
    direction_t::Vector{Float64},
)
    dcglp = oracle.dcglp

    if oracle.param.reuse_dcglp && haskey(dcglp, :conx) && haskey(dcglp, :cont) &&
       length(dcglp[:conx]) == length(direction_x) && length(dcglp[:cont]) == length(direction_t)
        tau = dcglp[:tau]
        for j in eachindex(direction_x)
            set_normalized_coefficient(dcglp[:conx][j], tau, direction_x[j])
            set_normalized_rhs(dcglp[:conx][j], x_value[j])
        end
        for j in eachindex(direction_t)
            set_normalized_coefficient(dcglp[:cont][j], tau, direction_t[j])
            set_normalized_rhs(dcglp[:cont][j], t_value[j])
        end
        return nothing
    end

    delete_registered_constraints!(dcglp, :conx)
    delete_registered_constraints!(dcglp, :cont)

    @constraint(
        dcglp,
        conx[j in eachindex(direction_x)],
        dcglp[:omega_x][1, j] + dcglp[:omega_x][2, j] + direction_x[j] * dcglp[:tau] == x_value[j],
    )
    @constraint(
        dcglp,
        cont[j in eachindex(direction_t)],
        dcglp[:omega_t][1, j] + dcglp[:omega_t][2, j] + direction_t[j] * dcglp[:tau] == t_value[j],
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
