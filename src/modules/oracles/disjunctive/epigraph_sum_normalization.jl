function add_normalization_constraint!(
    dcglp::Model,
    ::EpigraphSumNormalization,
    master::AbstractMaster,
    ::SplitOracleParam{EpigraphSumNormalization},
)
    @variable(dcglp, tau[1:master.dim_t])

    @objective(dcglp, Min, sum(tau))

    @constraint(dcglp, conx[j in 1:master.dim_x], dcglp[:omega_x][1, j] + dcglp[:omega_x][2, j] == 0.0)
    @constraint(dcglp, cont[j in 1:master.dim_t], dcglp[:omega_t][1, j] + dcglp[:omega_t][2, j] - tau[j] == 0.0)
end

function update_dcglp_for_candidate!(::EpigraphSumNormalization, oracle::SplitOracle{EpigraphSumNormalization}, x_value::Vector{Float64}, ::Vector{Float64})
    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
end

dcglp_tau_value(::EpigraphSumNormalization, dcglp::Model) = value.(dcglp[:tau])

function update_dcglp_upper_bound_and_gap!(
    ::EpigraphSumNormalization,
    state::DcglpState,
    log::DcglpLog,
    ::Vector{Float64},
    t_value::Vector{Float64},
)
    fill_dcglp_omega_t_estimates!(state, t_value)
    all(f_i -> !any(isnan, f_i), state.f_x) || return nothing
    update_upper_bound_and_gap!(state, log, (t1, t2) -> sum(t1) + sum(t2))
end

function has_dcglp_disjunctive_cut(::EpigraphSumNormalization, current_lb::Float64, t_value::Vector{Float64}, zero_tol::Float64)
    return current_lb >= sum(t_value) + zero_tol
end

function build_dcglp_disjunctive_cut(
    ::EpigraphSumNormalization,
    dcglp::Model,
    common::SplitOracleParam{EpigraphSumNormalization},
    ::Float64,
    ::Vector{Float64},
    t_value::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int},
)
    gamma_x = dual.(dcglp[:conx])
    gamma_0 = dual(dcglp[:con0])
    gamma_x, gamma_0 = apply_lift_or_strengthen(
        dcglp, gamma_x, zero_indices, one_indices;
        lift = common.lift, strengthen = common.strengthened,
        zero_tol = common.zero_tol, gamma_0 = gamma_0,
    )
    return Hyperplane(gamma_x, fill(-1.0, length(t_value)), gamma_0)
end
