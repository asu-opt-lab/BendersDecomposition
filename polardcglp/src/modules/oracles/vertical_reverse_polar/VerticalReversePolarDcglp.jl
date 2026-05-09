function BendersX.solve_dcglp!(
    oracle::VerticalReversePolarDCGLPOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int};
    start_time::Float64,
    time_limit::Float64,
    throw_typical_cuts_for_errors::Bool,
    include_disjunctive_cuts_to_hyperplanes::Bool,
)
    log = BendersX.DcglpLog()
    dcglp = oracle.dcglp
    typical_oracles = oracle.typical_oracles
    hyperplanes = Vector{BendersX.Hyperplane}()

    while true
        state = BendersX.DcglpState()
        state.total_time = @elapsed begin
            state.master_time = @elapsed begin
                JuMP.set_time_limit_sec(dcglp, BendersX.get_sec_remaining(log.start_time, time_limit))
                try
                    optimize!(dcglp)
                catch err
                    return fallback_typical_or_throw(
                        oracle,
                        x_value,
                        t_value,
                        start_time,
                        time_limit,
                        "VerticalReversePolarDCGLP master: unexpected error encountered when optimizing dcglp master: $(err)";
                        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
                    )
                end

                if is_solved_and_feasible(dcglp; allow_local = false, dual = true)
                    for i in 1:2
                        state.values[:ω_x][i] = value.(dcglp[:omega_x][i, :])
                        state.values[:ω_t][i] = value.(dcglp[:omega_t][i, :])
                        state.values[:ω_0][i] = value(dcglp[:omega_0][i])
                    end
                    state.values[:tau] = value(dcglp[:tau])
                    state.values[:sx] = Float64[]
                    state.LB = objective_value(dcglp)
                elseif termination_status(dcglp) == ALMOST_INFEASIBLE
                    return fallback_typical_or_throw(
                        oracle,
                        x_value,
                        t_value,
                        start_time,
                        time_limit,
                        "VerticalReversePolarDCGLP master: unexpected dcglp master termination status: $(termination_status(dcglp)); the problem is infeasible or dcglp encountered numerical issue";
                        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
                    )
                elseif termination_status(dcglp) == TIME_LIMIT
                    throw(BendersX.TimeLimitException("Time limit reached during VerticalReversePolarDCGLP solving"))
                else
                    throw(BendersX.UnexpectedModelStatusException("VerticalReversePolarDCGLP master: $(termination_status(dcglp))"))
                end
            end

            benders_cuts = Dict(1 => Vector{AffExpr}(), 2 => Vector{AffExpr}())
            for i in 1:2
                state.oracle_times[i] = @elapsed begin
                    if state.values[:ω_0][i] >= oracle.param.zero_tol
                        state.is_in_L[i], hyperplanes_i, state.f_x[i] = BendersX.generate_cuts(
                            typical_oracles[i],
                            clamp.(state.values[:ω_x][i] ./ state.values[:ω_0][i], 0.0, 1.0),
                            state.values[:ω_t][i] ./ state.values[:ω_0][i];
                            tol_normalize = state.values[:ω_0][i],
                            time_limit = BendersX.get_sec_remaining(log.start_time, time_limit),
                        )

                        if !state.is_in_L[i]
                            for k in 1:2
                                append!(
                                    benders_cuts[i],
                                    BendersX.hyperplanes_to_expression(
                                        dcglp,
                                        hyperplanes_i,
                                        dcglp[:omega_x][k, :],
                                        dcglp[:omega_t][k, :],
                                        dcglp[:omega_0][k],
                                    ),
                                )
                            end
                            if oracle.param.add_benders_cuts_to_master != 0
                                add_violated = oracle.param.add_benders_cuts_to_master == 2
                                append!(
                                    hyperplanes,
                                    BendersX.select_top_fraction(
                                        hyperplanes_i,
                                        h -> BendersX.evaluate_violation(h, x_value, t_value; zero_tol = oracle.param.zero_tol),
                                        oracle.param.fraction_of_benders_cuts_to_master;
                                        add_only_violated_cuts = add_violated,
                                    ),
                                )
                            end
                        end
                    else
                        state.is_in_L[i] = true
                        state.f_x[i] = zeros(length(t_value))
                    end
                end
            end

            for i in 1:2
                state.omega_t_[i] =
                    state.is_in_L[i] ? state.values[:ω_t][i] :
                    any(isnan, state.f_x[i]) ? fill(NaN, length(t_value)) :
                    state.f_x[i] * state.values[:ω_0][i]
            end
            if all(f_i -> !any(isnan, f_i), state.f_x)
                BendersX.update_upper_bound_and_gap!(state, log, (t1, t2) -> maximum(t1 .+ t2))
            end

            BendersX.record_iteration!(log, state)
        end

        oracle.param.dcglp_param.verbose && print_vertical_reverse_polar_iteration_info(state, log)
        BendersX.check_lb_improvement!(state, log; zero_tol = oracle.param.zero_tol)
        BendersX.is_terminated(state, log, oracle.param.dcglp_param, time_limit) && break

        BendersX.add_constraints(dcglp, :con_benders, [benders_cuts[1]; benders_cuts[2]])
    end

    current_lb = log.iterations[end].LB
    if current_lb >= oracle.param.zero_tol
        cut = generate_vertical_reverse_polar_disjunctive_cut(
            dcglp,
            current_lb,
            x_value,
            t_value,
            length(t_value),
            zero_indices,
            one_indices;
            strengthen = oracle.param.strengthened,
            lift = oracle.param.lift,
            zero_tol = oracle.param.zero_tol,
        )
        oracle.param.dcglp_param.verbose &&
            print_disjunctive_cut(oracle, cut, x_value, t_value; zero_tol = oracle.param.zero_tol)
        store_dcglp_disjunctive_cut!(
            oracle,
            cut,
            hyperplanes,
            include_disjunctive_cuts_to_hyperplanes,
        )
        return false, hyperplanes, fill(Inf, length(t_value))
    end

    return BendersX.generate_cuts(
        oracle.typical_oracles[1],
        x_value,
        t_value;
        time_limit = BendersX.get_sec_remaining(start_time, time_limit),
    )
end

function generate_vertical_reverse_polar_disjunctive_cut(
    dcglp::Model,
    current_lb::Float64,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    dim_t::Int,
    zero_indices::Vector{Int},
    one_indices::Vector{Int};
    strengthen::Bool = false,
    lift::Bool = false,
    zero_tol::Float64 = 1e-9,
)
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    gamma_0 = current_lb - dot(gamma_x, x_value) - dot(gamma_t, t_value)

    if lift && (!isempty(zero_indices) || !isempty(one_indices))
        zeta_k = !isempty(zero_indices) ? dual.(dcglp[:con_zeta][1, :]) : Float64[]
        zeta_v = !isempty(zero_indices) ? dual.(dcglp[:con_zeta][2, :]) : Float64[]
        xi_k = !isempty(one_indices) ? dual.(dcglp[:con_xi][1, :]) : Float64[]
        xi_v = !isempty(one_indices) ? dual.(dcglp[:con_xi][2, :]) : Float64[]

        gamma_0 -= sum(max.(xi_k, xi_v))
        lifted_gamma_x = -gamma_x
        lifted_gamma_x[zero_indices] .= -gamma_x[zero_indices] .+ max.(zeta_k, zeta_v)
        lifted_gamma_x[one_indices] .= -gamma_x[one_indices] .- max.(xi_k, xi_v)

        if strengthen
            sigma = Dict(1 => dual(dcglp[:con_split_kappa]), 2 => dual(dcglp[:con_split_nu]))
            delta_1 = dual.(dcglp[:condelta][1, :])
            delta_2 = dual.(dcglp[:condelta][2, :])
            delta_1[zero_indices] .+= -zeta_k .+ max.(zeta_k, zeta_v)
            delta_2[zero_indices] .+= -zeta_v .+ max.(zeta_k, zeta_v)
            delta = Dict(1 => delta_1, 2 => delta_2)
            lifted_gamma_x = BendersX.strengthening!(lifted_gamma_x, sigma, delta; zero_tol = zero_tol)
        end

        gamma_x = -lifted_gamma_x
    elseif strengthen
        sigma = Dict(1 => dual(dcglp[:con_split_kappa]), 2 => dual(dcglp[:con_split_nu]))
        delta = Dict(1 => dual.(dcglp[:condelta][1, :]), 2 => dual.(dcglp[:condelta][2, :]))
        gamma_x = -BendersX.strengthening!(-gamma_x, sigma, delta; zero_tol = zero_tol)
    end

    length(gamma_t) == dim_t ||
        throw(DimensionMismatch("`cont` dual has length $(length(gamma_t)) but expected $(dim_t)."))

    return BendersX.Hyperplane(gamma_x, gamma_t, gamma_0)
end

function fallback_typical_or_throw(
    oracle::VerticalReversePolarDCGLPOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
    msg::String;
    throw_typical_cuts_for_errors::Bool,
)
    if throw_typical_cuts_for_errors
        @warn msg
        return BendersX.generate_cuts(oracle.typical_oracles[1], x_value, t_value; time_limit = BendersX.get_sec_remaining(start_time, time_limit))
    end
    throw(BendersX.UnexpectedModelStatusException(msg))
end

function print_vertical_reverse_polar_iteration_info(state::BendersX.DcglpState, log::BendersX.DcglpLog)
    @printf(
        "   Iter: %4d | LB: %8.4f | UB: %8.4f | Gap: %6.2f%% | UB_k: %8.2f | UB_v: %8.2f | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
        log.n_iter,
        state.LB,
        state.UB,
        state.gap,
        maximum(state.omega_t_[1]),
        maximum(state.omega_t_[2]),
        state.master_time,
        state.oracle_times[1],
        state.oracle_times[2],
    )
    return nothing
end

function print_disjunctive_cut(
    oracle::VerticalReversePolarDCGLPOracle,
    cut::BendersX.Hyperplane,
    x_value::Vector{Float64},
    t_value::Vector{Float64};
    zero_tol::Float64 = 1e-10,
)
    println("   Disjunctive cut:")
    println("      " * format_hyperplane(cut; zero_tol = zero_tol))
    if isa(oracle.param.split_index_selection_rule, BendersX.SimpleSplit)
        @printf("   Split index: x[%d]\n", get_split_index(oracle))
    end
    @printf("   Cut violation at current point: %.6f\n", hyperplane_violation(cut, x_value, t_value))
    return nothing
end
