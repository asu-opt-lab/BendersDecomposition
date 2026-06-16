function add_normalization_constraint!(
    ::VerticalReversePolarNormalization,
    dcglp::Model,
    tau::VariableRef,
    sx::AbstractVector{VariableRef},
    st::AbstractVector{VariableRef},
)
    @constraint(dcglp, con_reverse_polar_x[j in eachindex(sx)], sx[j] + 0.0 * tau == 0.0)
    @constraint(dcglp, con_reverse_polar_t[j in eachindex(st)], st[j] + 0.0 * tau == 0.0)
end

function update_dcglp_for_candidate!(::VerticalReversePolarNormalization, oracle::SplitOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    set_normalized_rhs.(oracle.dcglp[:cont], t_value)
    update_reverse_polar_constraints!(oracle, zeros(length(x_value)), fill(-1.0, length(t_value)))
end

function update_dcglp_upper_bound_and_gap!(
    ::Union{VerticalReversePolarNormalization, DirectionalReversePolarNormalization},
    state::DcglpState,
    log::DcglpLog,
    ::Vector{Float64},
    t_value::Vector{Float64},
)
    fill_dcglp_omega_t_estimates!(state, t_value)
    previous_ub = log.n_iter >= 1 ? log.iterations[end].UB : Inf
    state.UB = state.is_in_L[1] && state.is_in_L[2] ? min(previous_ub, state.LB) : previous_ub
    state.gap =
        isfinite(state.UB) ?
        max(state.UB - state.LB, 0.0) / max(abs(state.UB), 1.0) * 100 :
        Inf
end
