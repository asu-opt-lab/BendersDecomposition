function BendersX.solve_dcglp!(
    oracle::DirectionalPolarDCGLPOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    direction_x::Vector{Float64},
    direction_t::Vector{Float64},
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
        state.values[:oracle_statuses] = fill(:inactive, 2)
        state.values[:oracle_gaps] = fill(0.0, 2)
        state.values[:oracle_scaled_gaps] = fill(0.0, 2)

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
                        "DirectionalPolarDCGLP master: unexpected error encountered when optimizing dcglp master: $(err)";
                        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
                    )
                end

                if is_solved_and_feasible(dcglp; allow_local = false, dual = true)
                    for i in 1:2
                        state.values[:ω_x][i] = value.(dcglp[:omega_x][i, :])
                        state.values[:ω_t][i] = value.(dcglp[:omega_t][i, :])
                        state.values[:ω_0][i] = value(dcglp[:omega_0][i])
                    end
                    state.values[:tau] = value(dcglp[:tau_var])
                    state.values[:sx] = Float64[]
                    state.LB = objective_value(dcglp)
                elseif termination_status(dcglp) == ALMOST_INFEASIBLE
                    return fallback_typical_or_throw(
                        oracle,
                        x_value,
                        t_value,
                        start_time,
                        time_limit,
                        "DirectionalPolarDCGLP master: unexpected dcglp master termination status: $(termination_status(dcglp)); the problem is infeasible or dcglp encountered numerical issue";
                        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
                    )
                elseif termination_status(dcglp) == TIME_LIMIT
                    throw(BendersX.TimeLimitException("Time limit reached during DirectionalPolarDCGLP solving"))
                else
                    throw(BendersX.UnexpectedModelStatusException("DirectionalPolarDCGLP master: $(termination_status(dcglp))"))
                end
            end

            benders_cuts = Dict(1 => Vector{AffExpr}(), 2 => Vector{AffExpr}())
            for i in 1:2
                state.oracle_times[i] = @elapsed begin
                    if state.values[:ω_0][i] >= oracle.param.zero_tol
                        t_block = state.values[:ω_t][i] ./ state.values[:ω_0][i]
                        state.is_in_L[i], hyperplanes_i, state.f_x[i] = BendersX.generate_cuts(
                            typical_oracles[i],
                            clamp.(state.values[:ω_x][i] ./ state.values[:ω_0][i], 0.0, 1.0),
                            t_block;
                            tol_normalize = state.values[:ω_0][i],
                            time_limit = BendersX.get_sec_remaining(log.start_time, time_limit),
                        )

                        state.values[:oracle_statuses][i] = state.is_in_L[i] ? :feasible : :infeasible
                        state.values[:oracle_gaps][i], state.values[:oracle_scaled_gaps][i] = compute_directional_oracle_gap(
                            state.is_in_L[i],
                            t_block,
                            state.f_x[i],
                            state.values[:ω_0][i],
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
                                        h -> hyperplane_violation(h, x_value, t_value),
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
            update_directional_upper_bound_and_gap!(state, log)
            BendersX.record_iteration!(log, state)
        end

        oracle.param.dcglp_param.verbose && print_directional_iteration_info(state, log)
        BendersX.check_lb_improvement!(state, log; zero_tol = oracle.param.zero_tol)
        BendersX.is_terminated(state, log, oracle.param.dcglp_param, time_limit) && break

        BendersX.add_constraints(dcglp, :con_benders, [benders_cuts[1]; benders_cuts[2]])
    end

    current_lb = log.iterations[end].LB
    if current_lb >= oracle.param.zero_tol
        cut = generate_directional_polar_cut(
            dcglp,
            x_value,
            t_value,
            direction_x,
            direction_t,
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

function generate_directional_polar_cut(
    dcglp::Model,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    direction_x::Vector{Float64},
    direction_t::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int};
    strengthen::Bool = false,
    lift::Bool = false,
    zero_tol::Float64 = 1e-9,
)
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    normalize_directional_duals!(gamma_x, gamma_t, direction_x, direction_t; zero_tol = zero_tol)

    gamma_0 = -dot(gamma_x, x_value) - dot(gamma_t, t_value) + value(dcglp[:tau_var])

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

    return BendersX.Hyperplane(gamma_x, gamma_t, gamma_0)
end

function normalize_directional_duals!(
    gamma_x::AbstractVector{Float64},
    gamma_t::AbstractVector{Float64},
    direction_x::Vector{Float64},
    direction_t::Vector{Float64};
    zero_tol::Float64,
)
    direction_value = dot(gamma_x, direction_x) + dot(gamma_t, direction_t)
    abs(direction_value) > zero_tol ||
        throw(BendersX.AlgorithmException("DirectionalPolarDCGLP cut normalization failed because the directional support is numerically zero."))
    gamma_x ./= direction_value
    gamma_t ./= direction_value
    return nothing
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

function fallback_typical_or_throw(
    oracle::DirectionalPolarDCGLPOracle,
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

function update_directional_upper_bound_and_gap!(state::BendersX.DcglpState, log::BendersX.DcglpLog)
    previous_ub = log.n_iter >= 1 ? log.iterations[end].UB : 1.0
    state.UB = state.is_in_L[1] && state.is_in_L[2] ? min(previous_ub, state.LB) : previous_ub
    state.gap =
        isfinite(state.UB) ?
        max(state.UB - state.LB, 0.0) / max(abs(state.UB), 1.0) * 100 :
        Inf
    return nothing
end

function print_directional_iteration_info(state::BendersX.DcglpState, log::BendersX.DcglpLog)
    @printf(
        "   Iter: %4d | LB: %12.8f | UB: %12.8f | Gap: %8.4f%% | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
        log.n_iter,
        state.LB,
        state.UB,
        state.gap,
        state.master_time,
        state.oracle_times[1],
        state.oracle_times[2],
    )
    @printf(
        "      Oracle status | k: %-10s gap: %8.4f scaled: %8.4f | v: %-10s gap: %8.4f scaled: %8.4f\n",
        String(state.values[:oracle_statuses][1]),
        state.values[:oracle_gaps][1],
        state.values[:oracle_scaled_gaps][1],
        String(state.values[:oracle_statuses][2]),
        state.values[:oracle_gaps][2],
        state.values[:oracle_scaled_gaps][2],
    )
    return nothing
end

function print_disjunctive_cut(
    oracle::DirectionalPolarDCGLPOracle,
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
