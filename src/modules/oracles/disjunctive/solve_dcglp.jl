"""
    should_fallback_typical(normalization, oracle, x_value, t_value) -> Bool

If `true`, `generate_cuts` short-circuits and delegates to the first typical
oracle without solving the DCGLP. Default implementation returns `false`.
The directional normalization uses this hook both to skip the DCGLP when the
candidate point coincides with the core point and to cache the direction
vector for later cut extraction.
"""
should_fallback_typical(::AbstractDisjunctiveNormalization, ::SplitOracle, ::Vector{Float64}, ::Vector{Float64}) = false

function solve_dcglp_loop!(
    oracle::SplitOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int};
    start_time::Float64,
    time_limit::Float64,
    throw_typical_cuts_for_errors::Bool,
    include_disjunctive_cuts_to_hyperplanes::Bool,
)
    normalization = oracle.param.normalization
    normalization_name = string(nameof(typeof(normalization)))
    log = DcglpLog()
    log.start_time = start_time

    dcglp = oracle.dcglp
    hyperplanes = Hyperplane[]
    reference_t = update_dcglp_reference_t!(normalization, oracle, x_value, t_value, start_time, time_limit)

    while true
        state = initialize_dcglp_state(normalization)
        benders_cuts = Dict(1 => AffExpr[], 2 => AffExpr[])

        state.total_time = @elapsed begin
            state.master_time = @elapsed begin
                set_time_limit_sec(dcglp, get_sec_remaining(log.start_time, time_limit))
                try
                    optimize!(dcglp)
                catch err
                    return fallback_typical_or_throw(
                        oracle, x_value, t_value, start_time, time_limit,
                        "$(normalization_name) master: unexpected error encountered when optimizing dcglp master: $(err)";
                        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
                    )
                end

                if is_solved_and_feasible(dcglp; allow_local = false, dual = true)
                    read_dcglp_solution!(oracle, state)
                elseif termination_status(dcglp) == ALMOST_INFEASIBLE
                    return fallback_typical_or_throw(
                        oracle, x_value, t_value, start_time, time_limit,
                        "$(normalization_name) master: unexpected dcglp master termination status: $(termination_status(dcglp)); the problem is infeasible or dcglp encountered numerical issue";
                        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
                    )
                elseif termination_status(dcglp) == TIME_LIMIT
                    throw(TimeLimitException("Time limit reached during $(normalization_name) solving"))
                else
                    return fallback_typical_or_throw(
                        oracle, x_value, t_value, start_time, time_limit,
                        "$(normalization_name) master: termination status is $(termination_status(dcglp))";
                        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
                    )
                end
            end

            collect_dcglp_benders_cuts!(oracle, state, benders_cuts, hyperplanes, x_value, t_value, log, time_limit)
            update_dcglp_upper_bound_and_gap!(normalization, state, log, reference_t, t_value)
            record_iteration!(log, state)
        end

        oracle.param.dcglp_param.verbose && print_dcglp_iteration_info(normalization, state, log)
        check_lb_improvement!(state, log; zero_tol = oracle.param.zero_tol)
        is_terminated(state, log, oracle.param.dcglp_param, time_limit) && break

        cuts_to_add = [benders_cuts[1]; benders_cuts[2]]
        !isempty(cuts_to_add) && add_constraints(dcglp, :con_benders, cuts_to_add)
    end

    current_lb = log.iterations[end].LB
    if current_lb >= oracle.param.zero_tol
        cut = build_dcglp_disjunctive_cut(normalization, dcglp, oracle.param, current_lb, x_value, t_value, zero_indices, one_indices)
        oracle.param.dcglp_param.verbose && print_disjunctive_cut(oracle, cut, x_value, t_value; zero_tol = oracle.param.zero_tol)
        store_dcglp_disjunctive_cut!(oracle, cut, hyperplanes, include_disjunctive_cuts_to_hyperplanes)
        return false, hyperplanes, fill(Inf, length(t_value))
    end

    return generate_cuts(
        oracle.typical_oracles[1],
        x_value,
        t_value;
        time_limit = get_sec_remaining(start_time, time_limit),
    )
end

function read_dcglp_solution!(oracle::SplitOracle, state::DcglpState)
    dcglp = oracle.dcglp
    for i in 1:2
        state.values[:ω_x][i] = value.(dcglp[:omega_x][i, :])
        state.values[:ω_t][i] = value.(dcglp[:omega_t][i, :])
        state.values[:ω_0][i] = value(dcglp[:omega_0][i])
    end
    tau_value = value(dcglp[:tau])
    state.values[:tau] = tau_value
    state.values[:sx] = value.(dcglp[:sx])
    state.LB = tau_value
end

function collect_dcglp_benders_cuts!(
    oracle::SplitOracle,
    state::DcglpState,
    benders_cuts::Dict{Int, Vector{AffExpr}},
    hyperplanes::Vector{Hyperplane},
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    log::DcglpLog,
    time_limit::Float64,
)
    dcglp = oracle.dcglp
    normalization = oracle.param.normalization
    for i in 1:2
        state.oracle_times[i] = @elapsed begin
            if state.values[:ω_0][i] >= oracle.param.zero_tol
                t_block = state.values[:ω_t][i] ./ state.values[:ω_0][i]
                state.is_in_L[i], hyperplanes_i, state.f_x[i] = generate_cuts(
                    oracle.typical_oracles[i],
                    clamp.(state.values[:ω_x][i] ./ state.values[:ω_0][i], 0.0, 1.0),
                    t_block;
                    tol_normalize = state.values[:ω_0][i],
                    time_limit = get_sec_remaining(log.start_time, time_limit),
                )
                record_dcglp_oracle_result!(normalization, state, i, t_block)

                if !state.is_in_L[i]
                    for k in 1:2
                        append!(
                            benders_cuts[i],
                            hyperplanes_to_expression(
                                dcglp,
                                hyperplanes_i,
                                dcglp[:omega_x][k, :],
                                dcglp[:omega_t][k, :],
                                dcglp[:omega_0][k],
                            ),
                        )
                    end
                    append_selected_benders_cuts_to_master!(
                        oracle, hyperplanes, hyperplanes_i, x_value, t_value,
                    )
                end
            else
                state.is_in_L[i] = true
                state.f_x[i] = zeros(length(t_value))
            end
        end
    end
end

function append_selected_benders_cuts_to_master!(
    oracle::SplitOracle,
    hyperplanes::Vector{Hyperplane},
    candidate_hyperplanes::Vector{Hyperplane},
    x_value::Vector{Float64},
    t_value::Vector{Float64},
)
    oracle.param.add_benders_cuts_to_master == 0 && return nothing

    add_violated = oracle.param.add_benders_cuts_to_master == 2
    append!(
        hyperplanes,
        select_top_fraction(
            candidate_hyperplanes,
            h -> evaluate_violation(h, x_value, t_value),
            oracle.param.fraction_of_benders_cuts_to_master;
            add_only_violated_cuts = add_violated,
        ),
    )
end

function fill_dcglp_omega_t_estimates!(state::DcglpState, t_value::Vector{Float64})
    for i in 1:2
        state.omega_t_[i] =
            state.is_in_L[i] ? state.values[:ω_t][i] :
            any(isnan, state.f_x[i]) ? fill(NaN, length(t_value)) :
            state.f_x[i] * state.values[:ω_0][i]
    end
end
function build_dcglp_disjunctive_cut(
    ::AbstractDisjunctiveNormalization,
    dcglp::Model,
    common::SplitOracleParam,
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
    return Hyperplane(gamma_x, gamma_t, gamma_0)
end
