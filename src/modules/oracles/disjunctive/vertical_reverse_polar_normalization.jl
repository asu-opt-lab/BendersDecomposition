function add_normalization_constraint!(
    dcglp::Model,
    ::VerticalReversePolarNormalization,
    master::AbstractMaster,
    ::SplitOracleParam{VerticalReversePolarNormalization},
)
    @variable(dcglp, tau)

    @objective(dcglp, Min, tau)

    @constraint(dcglp, conx[j in 1:master.dim_x], dcglp[:omega_x][1, j] + dcglp[:omega_x][2, j] == 0.0)
    @constraint(dcglp, cont[j in 1:master.dim_t], dcglp[:omega_t][1, j] + dcglp[:omega_t][2, j] - tau == 0.0)
end

function update_dcglp_for_candidate!(::VerticalReversePolarNormalization, oracle::SplitOracle{VerticalReversePolarNormalization}, x_value::Vector{Float64}, t_value::Vector{Float64})
    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    set_normalized_rhs.(oracle.dcglp[:cont], t_value)
end

dcglp_tau_value(::VerticalReversePolarNormalization, dcglp::Model) = value(dcglp[:tau])

function update_dcglp_upper_bound_and_gap!(
    ::VerticalReversePolarNormalization,
    state::DcglpState,
    log::DcglpLog,
    ::Vector{Float64},
    t_value::Vector{Float64},
)
    fill_dcglp_omega_t_estimates!(state, t_value)
    all(f_i -> !any(isnan, f_i), state.f_x) || return nothing
    update_upper_bound_and_gap!(state, log, (t1, t2) -> maximum(t1 .+ t2))
end

function has_dcglp_disjunctive_cut(::VerticalReversePolarNormalization, current_lb::Float64, ::Vector{Float64}, zero_tol::Float64)
    return current_lb >= zero_tol
end

function build_dcglp_disjunctive_cut(
    ::VerticalReversePolarNormalization,
    dcglp::Model,
    common::SplitOracleParam{VerticalReversePolarNormalization},
    ::Float64,
    ::Vector{Float64},
    t_value::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int},
)
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    gamma_0 = dual(dcglp[:con0])
    gamma_x, gamma_0 = apply_lift_or_strengthen(
        dcglp, gamma_x, zero_indices, one_indices;
        lift = common.lift, strengthen = common.strengthened,
        zero_tol = common.zero_tol, gamma_0 = gamma_0,
    )
    length(gamma_t) == length(t_value) ||
        throw(DimensionMismatch("`cont` dual has length $(length(gamma_t)) but expected $(length(t_value))."))
    return Hyperplane(gamma_x, gamma_t, gamma_0)
end
