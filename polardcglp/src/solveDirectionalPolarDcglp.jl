function update_directional_no_improvement!(
    no_improvement::Int,
    current_lb::Float64,
    prev_lb::Float64;
    zero_tol::Float64 = 1e-8,
    tol_imprv::Float64 = 1e-4,
)
    lb_improvement =
        abs(prev_lb) < zero_tol ? abs(current_lb - prev_lb) :
        abs((current_lb - prev_lb) / prev_lb) * 100
    return lb_improvement < tol_imprv ? no_improvement + 1 : 0
end

function print_directional_iteration_info(
    iteration::Int,
    lb::Float64,
    ub::Float64,
    gap::Float64,
    master_time::Float64,
    oracle_times::Vector{Float64},
)
    @printf(
        "   Iter: %4d | LB: %8.4f | UB: %8.4f | Gap: %6.2f%% | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
        iteration,
        lb,
        ub,
        gap,
        master_time,
        oracle_times[1],
        oracle_times[2],
    )
    return nothing
end

function print_directional_stop_reason(msg::String)
    @printf("   Stop: %s\n", msg)
    return nothing
end

function format_sparse_terms(coeffs::SparseVector{Float64, Int}, var_name::String; zero_tol::Float64 = 1e-10)
    terms = String[]
    for p in sortperm(coeffs.nzind)
        idx = coeffs.nzind[p]
        val = coeffs.nzval[p]
        abs(val) <= zero_tol && continue
        sign = val >= 0 ? "+" : "-"
        push!(terms, @sprintf("%s %.6g %s[%d]", sign, abs(val), var_name, idx))
    end
    return terms
end

function format_hyperplane(h::BendersX.Hyperplane; zero_tol::Float64 = 1e-10)
    pieces = String[]
    if abs(h.a_0) > zero_tol || (nnz(h.a_x) == 0 && nnz(h.a_t) == 0)
        push!(pieces, @sprintf("%.6g", h.a_0))
    else
        push!(pieces, "0")
    end
    append!(pieces, format_sparse_terms(h.a_x, "x"; zero_tol = zero_tol))
    append!(pieces, format_sparse_terms(h.a_t, "t"; zero_tol = zero_tol))
    return join(pieces, " ") * " <= 0"
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

hyperplane_violation(h::BendersX.Hyperplane, x_value::Vector{Float64}, t_value::Vector{Float64}) =
    h.a_0 + dot(h.a_x, x_value) + dot(h.a_t, t_value)

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

function optimize_directional_polar_dcglp!(
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
    dcglp = oracle.dcglp
    master_hyperplanes = BendersX.Hyperplane[]

    # tau = 1 corresponds to the core point, which is assumed feasible.
    best_lb = 0.0
    best_ub = 1.0
    prev_lb = 0.0
    current_lb = 0.0
    no_improvement = 0
    iteration = 0

    while true
        iteration += 1
        JuMP.set_time_limit_sec(dcglp, BendersX.get_sec_remaining(start_time, time_limit))

        master_time = 0.0
        try
            master_time = @elapsed optimize!(dcglp)
        catch err
            return fallback_typical_or_throw(
                oracle,
                x_value,
                t_value,
                start_time,
                time_limit,
                "DirectionalPolarDCGLP master failed during optimization: $(err)";
                throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
            )
        end

        status = termination_status(dcglp)
        if status == TIME_LIMIT
            throw(BendersX.TimeLimitException("Time limit reached during DirectionalPolarDCGLP solving."))
        elseif status != OPTIMAL
            return fallback_typical_or_throw(
                oracle,
                x_value,
                t_value,
                start_time,
                time_limit,
                "DirectionalPolarDCGLP master terminated with unexpected status $(status).";
                throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
            )
        end

        current_lb = objective_value(dcglp)
        if iteration >= 2
            no_improvement = update_directional_no_improvement!(
                no_improvement,
                current_lb,
                prev_lb;
                zero_tol = oracle.param.zero_tol,
            )
        end
        best_lb = max(best_lb, current_lb)
        prev_lb = current_lb

        violated_cuts = AffExpr[]
        current_points_feasible = true
        oracle_times = zeros(2)

        for i in 1:2
            omega_0 = value(dcglp[:omega_0][i])
            if omega_0 < oracle.param.zero_tol
                continue
            end

            x_block = clamp.(value.(dcglp[:omega_x][i, :]) ./ omega_0, 0.0, 1.0)
            t_block = value.(dcglp[:omega_t][i, :]) ./ omega_0

            f_x_i = similar(t_value)
            is_in_L = false
            hyperplanes_i = BendersX.Hyperplane[]
            oracle_times[i] = @elapsed begin
                is_in_L, hyperplanes_i, f_x_i = BendersX.generate_cuts(
                    oracle.typical_oracles[i],
                    x_block,
                    t_block;
                    tol_normalize = omega_0,
                    time_limit = BendersX.get_sec_remaining(start_time, time_limit),
                )
            end

            if !is_in_L
                current_points_feasible = false
                for block_idx in 1:2
                    append!(violated_cuts, BendersX.hyperplanes_to_expression(dcglp, hyperplanes_i, dcglp[:omega_x][block_idx, :], dcglp[:omega_t][block_idx, :], dcglp[:omega_0][block_idx]))
                end

                if oracle.param.add_benders_cuts_to_master != 0
                    add_only_violated = oracle.param.add_benders_cuts_to_master == 2
                    selected = BendersX.select_top_fraction(
                        hyperplanes_i,
                        h -> hyperplane_violation(h, x_value, t_value),
                        oracle.param.fraction_of_benders_cuts_to_master;
                        add_only_violated_cuts = add_only_violated,
                    )
                    append!(master_hyperplanes, selected)
                end
            end
        end

        if current_points_feasible
            best_ub = min(best_ub, current_lb)
        end

        gap =
            isfinite(best_ub) ?
            max(best_ub - current_lb, 0.0) / max(abs(best_ub), 1.0) * 100 :
            Inf

        oracle.param.dcglp_param.verbose &&
            print_directional_iteration_info(iteration, current_lb, best_ub, gap, master_time, oracle_times)

        if current_points_feasible
            if current_lb > oracle.param.zero_tol
                oracle.param.dcglp_param.verbose &&
                    print_directional_stop_reason("directional boundary found before the current point; generating a disjunctive cut")
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
                append_current_disjunctive_cut!(oracle, cut)
                if include_disjunctive_cuts_to_hyperplanes
                    push!(master_hyperplanes, cut)
                end
                return false, master_hyperplanes, fill(Inf, length(t_value))
            end

            oracle.param.dcglp_param.verbose &&
                print_directional_stop_reason("current point remains on the feasible segment; falling back to typical cuts")
            return BendersX.generate_cuts(oracle.typical_oracles[1], x_value, t_value; time_limit = BendersX.get_sec_remaining(start_time, time_limit))
        end

        if gap <= oracle.param.dcglp_param.gap_tolerance
            oracle.param.dcglp_param.verbose &&
                print_directional_stop_reason("dcglp gap tolerance reached")
            break
        end

        if no_improvement >= oracle.param.dcglp_param.halt_limit
            oracle.param.dcglp_param.verbose &&
                print_directional_stop_reason("halt limit reached due to insufficient LB improvement")
            break
        end

        if iteration >= oracle.param.dcglp_param.iter_limit
            oracle.param.dcglp_param.verbose &&
                print_directional_stop_reason("iteration limit reached")
            break
        end

        elapsed_time = time() - start_time
        if elapsed_time >= time_limit || elapsed_time >= oracle.param.dcglp_param.time_limit
            oracle.param.dcglp_param.verbose &&
                print_directional_stop_reason("time limit reached")
            break
        end

        BendersX.add_constraints(dcglp, :con_benders, violated_cuts)
    end

    if current_lb > oracle.param.zero_tol
        oracle.param.dcglp_param.verbose &&
            print_directional_stop_reason("LB sufficient; generating a disjunctive cut")
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
        append_current_disjunctive_cut!(oracle, cut)
        if include_disjunctive_cuts_to_hyperplanes
            push!(master_hyperplanes, cut)
        end
        return false, master_hyperplanes, fill(Inf, length(t_value))
    end

    return fallback_typical_or_throw(
        oracle,
        x_value,
        t_value,
        start_time,
        time_limit,
        "DirectionalPolarDCGLP stopped before certifying the directional boundary; returning typical cuts.";
        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
    )
end
