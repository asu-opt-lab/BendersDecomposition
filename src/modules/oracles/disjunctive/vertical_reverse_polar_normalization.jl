function build_dcglp(
    master::AbstractMaster,
    param::SplitOracleParam,
    ::VerticalReversePolarNormalization,
)
    dcglp = Model(param.dcglp_param.optimizer)
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

    @variable(dcglp, tau)

    @objective(dcglp, Min, tau)

    @constraint(dcglp, conx[j in 1:master.dim_x], dcglp[:omega_x][1, j] + dcglp[:omega_x][2, j] == 0.0)
    @constraint(dcglp, cont[j in 1:master.dim_t], dcglp[:omega_t][1, j] + dcglp[:omega_t][2, j] - tau == 0.0)

    return dcglp
end

function update_dcglp_for_candidate!(::VerticalReversePolarNormalization, oracle::SplitOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
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

