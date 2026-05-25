function build_strategy_dcglp(::SimplexNormStrategy, master::AbstractMaster, param::DcglpOracleParam{SimplexNormStrategy})
    dcglp = Model(param.dcglp_param.optimizer)

    @variable(dcglp, tau[1:master.dim_t])

    @objective(dcglp, Min, sum(tau))

    build_dcglp_skeleton!(dcglp, master; omega_0_nonneg = true)

    @constraint(dcglp, conx[j in 1:master.dim_x], dcglp[:omega_x][1, j] + dcglp[:omega_x][2, j] == 0.0)
    @constraint(dcglp, cont[j in 1:master.dim_t], dcglp[:omega_t][1, j] + dcglp[:omega_t][2, j] - tau[j] == 0.0)

    return dcglp
end

function prepare_dcglp_call!(::SimplexNormStrategy, oracle::SimplexNormOracle, x_value::Vector{Float64}, ::Vector{Float64})
    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    return :proceed
end

dcglp_tau_value(::SimplexNormStrategy, dcglp::Model) = value.(dcglp[:tau])

function update_dcglp_upper_bound_and_gap!(
    ::SimplexNormStrategy,
    state::DcglpState,
    log::DcglpLog,
    ::Vector{Float64},
    t_value::Vector{Float64},
)
    fill_dcglp_omega_t_estimates!(state, t_value)
    all(f_i -> !any(isnan, f_i), state.f_x) || return nothing
    update_upper_bound_and_gap!(state, log, (t1, t2) -> sum(t1) + sum(t2))
    return nothing
end

function has_dcglp_disjunctive_cut(::SimplexNormStrategy, current_lb::Float64, t_value::Vector{Float64}, zero_tol::Float64)
    return current_lb > sum(t_value) + zero_tol
end

function build_dcglp_disjunctive_cut(
    ::SimplexNormStrategy,
    dcglp::Model,
    common::DcglpOracleParam{SimplexNormStrategy},
    current_lb::Float64,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int},
)
    gamma_x = dual.(dcglp[:conx])
    gamma_0 = current_lb - dot(gamma_x, x_value)
    gamma_x, gamma_0 = apply_lift_or_strengthen(
        dcglp, gamma_x, zero_indices, one_indices;
        lift = common.lift, strengthen = common.strengthened,
        zero_tol = common.zero_tol, gamma_0 = gamma_0,
    )
    return Hyperplane(gamma_x, fill(-1.0, length(t_value)), gamma_0)
end

function print_dcglp_iteration_info(::SimplexNormStrategy, state::DcglpState, log::DcglpLog)
    @printf(
        "   Iter: %4d | LB: %8.4f | UB: %8.4f | Gap: %6.2f%% | UB_k: %8.2f | UB_v: %8.2f | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
        log.n_iter,
        state.LB,
        state.UB,
        state.gap,
        sum(state.omega_t_[1]),
        sum(state.omega_t_[2]),
        state.master_time,
        state.oracle_times[1],
        state.oracle_times[2],
    )
    return nothing
end
