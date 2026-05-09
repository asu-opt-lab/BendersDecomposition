function update_simplex_norm_old_no_improvement!(
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

function print_simplex_norm_old_iteration_info(
    iteration::Int,
    lb::Float64,
    ub::Float64,
    gap::Float64,
    ub_k::Float64,
    ub_v::Float64,
    master_time::Float64,
    oracle_times::Vector{Float64},
)
    @printf(
        "   Iter: %4d | LB: %8.4f | UB: %8.4f | Gap: %6.2f%% | UB_k: %8.2f | UB_v: %8.2f | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
        iteration,
        lb,
        ub,
        gap,
        ub_k,
        ub_v,
        master_time,
        oracle_times[1],
        oracle_times[2],
    )
    return nothing
end

function print_simplex_norm_old_stop_reason(msg::String)
    @printf("   Stop: %s\n", msg)
    return nothing
end

function print_disjunctive_cut(
    oracle::SimplexNormOldDCGLPOracle,
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

simplex_norm_old_master_t_value(t_value::Vector{Float64}) = sum(t_value)

function fallback_typical_or_throw(
    oracle::SimplexNormOldDCGLPOracle,
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

function print_dcglp_iis(
    oracle::SimplexNormOldDCGLPOracle,
    x_value::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int},
)
    dcglp = oracle.dcglp

    println("========== SimplexNormOldDCGLP IIS Diagnostic ==========")
    if !isempty(oracle.splits)
        phi, phi_0 = oracle.splits[end]
        j_split = findfirst(v -> v > 0.5, phi)
        if j_split !== nothing
            @printf("  Split index: j_split = %d, phi_0 = %.6g\n", j_split, phi_0)
            @printf("  x_value[j_split] = %.6g\n", x_value[j_split])
        end
    end
    if !isempty(zero_indices)
        println("  zero_indices (lifting): ", zero_indices)
    end
    if !isempty(one_indices)
        println("  one_indices (lifting): ", one_indices)
    end

    try
        JuMP.compute_conflict!(dcglp)
    catch err
        @warn "compute_conflict! failed: $(err)"
        println("======================================================")
        return
    end

    conflict_status = try
        MOI.get(dcglp, MOI.ConflictStatus())
    catch err
        @warn "Querying ConflictStatus failed: $(err)"
        println("======================================================")
        return
    end

    if conflict_status != MOI.CONFLICT_FOUND
        println("  Conflict refiner result: $(conflict_status) (no conflict reported)")
        println("======================================================")
        return
    end

    println("  Constraints in IIS:")
    n_in_conflict = 0
    for (F, S) in JuMP.list_of_constraint_types(dcglp)
        for con in JuMP.all_constraints(dcglp, F, S)
            cs = try
                MOI.get(dcglp, MOI.ConstraintConflictStatus(), con)
            catch
                continue
            end
            if cs == MOI.IN_CONFLICT
                n_in_conflict += 1
                nm = JuMP.name(con)
                label = isempty(nm) ? string(con) : nm
                println("    [$(F), $(S)] ", label)
            end
        end
    end
    @printf("  Total constraints in conflict: %d\n", n_in_conflict)
    println("======================================================")
    return nothing
end

function build_split_bound_cut(split_index::Int, dim_x::Int, dim_t::Int, side::Symbol)
    cut = BendersX.Hyperplane(dim_x, dim_t)
    if side == :x_ge_1
        cut.a_0 = 1.0
        cut.a_x[split_index] = -1.0
    elseif side == :x_le_0
        cut.a_x[split_index] = 1.0
    else
        throw(ArgumentError("Unsupported split bound cut side: $(side)"))
    end
    return cut
end

function restore_simplex_norm_old_conx!(dcglp::Model, x_value::Vector{Float64})
    delete_registered_constraints!(dcglp, :conx)
    @constraint(dcglp, conx[j in eachindex(x_value)], dcglp[:omega_x][1, j] + dcglp[:omega_x][2, j] == 0.0)
    JuMP.set_normalized_rhs.(dcglp[:conx], x_value)
    return nothing
end

function diagnose_infeasible_bound_cut(
    oracle::SimplexNormOldDCGLPOracle,
    x_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
)
    isa(oracle.param.split_index_selection_rule, BendersX.SimpleSplit) || return nothing
    isempty(oracle.splits) && return nothing

    dcglp = oracle.dcglp
    split_index = try
        get_split_index(oracle)
    catch
        return nothing
    end

    original_sense = JuMP.objective_sense(dcglp)
    original_objective = JuMP.objective_function(dcglp)
    omega_max = Dict{Int, Float64}()

    try
        delete_registered_constraints!(dcglp, :conx)

        for block_idx in 1:2
            JuMP.set_objective_sense(dcglp, MOI.MAX_SENSE)
            JuMP.set_objective_function(dcglp, dcglp[:omega_0][block_idx])
            JuMP.set_time_limit_sec(dcglp, BendersX.get_sec_remaining(start_time, time_limit))
            optimize!(dcglp)

            status = termination_status(dcglp)
            if status == OPTIMAL
                omega_max[block_idx] = value(dcglp[:omega_0][block_idx])
                @printf(
                    "  diagnose_infeasible_bound_cut: max omega_0[%d] = %.6g\n",
                    block_idx,
                    omega_max[block_idx],
                )
            elseif status == MOI.INFEASIBLE
                omega_max[block_idx] = 0.0
                @printf(
                    "  diagnose_infeasible_bound_cut: max omega_0[%d] infeasible -> %.6g\n",
                    block_idx,
                    omega_max[block_idx],
                )
            else
                @warn "SimplexNormOldDCGLP infeasible diagnostic for omega_0[$(block_idx)] returned status $(status); skipping temporary bound cut."
                return nothing
            end
        end
    finally
        JuMP.set_objective_sense(dcglp, original_sense)
        JuMP.set_objective_function(dcglp, original_objective)
        restore_simplex_norm_old_conx!(dcglp, x_value)
    end

    zero_tol = oracle.param.zero_tol
    omega_1 = get(omega_max, 1, 0.0)
    omega_2 = get(omega_max, 2, 0.0)

    if omega_2 <= zero_tol && omega_1 > zero_tol
        cut = build_split_bound_cut(split_index, length(x_value), length(dcglp[:tau]), :x_ge_1)
        println("  diagnose_infeasible_bound_cut: selected cut")
        println("    " * format_hyperplane(cut; zero_tol = zero_tol))
        return cut
    elseif omega_1 <= zero_tol && omega_2 > zero_tol
        cut = build_split_bound_cut(split_index, length(x_value), length(dcglp[:tau]), :x_le_0)
        println("  diagnose_infeasible_bound_cut: selected cut")
        println("    " * format_hyperplane(cut; zero_tol = zero_tol))
        return cut
    end

    if omega_1 <= zero_tol && omega_2 <= zero_tol
        @warn "SimplexNormOldDCGLP infeasible diagnostic found both split sides unsupported after removing conx; keeping fallback path."
    end
    return nothing
end

function generate_simplex_norm_old_disjunctive_cut(
    dcglp::Model,
    current_lb::Float64,
    x_value::Vector{Float64},
    dim_t::Int,
    zero_indices::Vector{Int},
    one_indices::Vector{Int};
    strengthen::Bool = false,
    lift::Bool = false,
    zero_tol::Float64 = 1e-9,
)
    gamma_x = dual.(dcglp[:conx])
    gamma_0 = current_lb - dot(gamma_x, x_value)

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

    return BendersX.Hyperplane(gamma_x, fill(-1.0, dim_t), gamma_0)
end

function optimize_simplex_norm_old_dcglp!(
    oracle::SimplexNormOldDCGLPOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int};
    start_time::Float64,
    time_limit::Float64,
    throw_typical_cuts_for_errors::Bool,
    include_disjunctive_cuts_to_hyperplanes::Bool,
)
    dcglp = oracle.dcglp
    master_hyperplanes = BendersX.Hyperplane[]

    best_lb = -Inf
    best_ub = Inf
    current_lb = -Inf
    no_improvement = 0
    iteration = 0
    prev_lb = -Inf

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
                "SimplexNormOldDCGLP master failed during optimization: $(err)";
                throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
            )
        end

        status = termination_status(dcglp)
        if status == TIME_LIMIT
            throw(BendersX.TimeLimitException("Time limit reached during SimplexNormOldDCGLP solving."))
        elseif status != OPTIMAL
            if status == MOI.INFEASIBLE
                print_dcglp_iis(oracle, x_value, zero_indices, one_indices)
                bound_cut = diagnose_infeasible_bound_cut(oracle, x_value, start_time, time_limit)
                if !isnothing(bound_cut)
                    oracle.param.dcglp_param.verbose &&
                        print_simplex_norm_old_stop_reason("dcglp infeasible; adding temporary split bound cut")
                    append_current_disjunctive_cut!(oracle, bound_cut)
                    if include_disjunctive_cuts_to_hyperplanes
                        push!(master_hyperplanes, bound_cut)
                    end
                    return false, master_hyperplanes, fill(Inf, length(t_value))
                end
            end
            return fallback_typical_or_throw(
                oracle,
                x_value,
                t_value,
                start_time,
                time_limit,
                "SimplexNormOldDCGLP master terminated with unexpected status $(status).";
                throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
            )
        end

        current_lb = objective_value(dcglp)
        if iteration >= 2
            no_improvement = update_simplex_norm_old_no_improvement!(
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
        omega_t_eval = [zeros(length(t_value)), zeros(length(t_value))]
        for i in 1:2
            omega_0 = value(dcglp[:omega_0][i])
            omega_t_value = value.(dcglp[:omega_t][i, :])
            omega_t_eval[i] .= omega_t_value
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

            if !any(isnan, f_x_i)
                omega_t_eval[i] .= is_in_L ? omega_t_value : f_x_i .* omega_0
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

        ub_k = sum(omega_t_eval[1])
        ub_v = sum(omega_t_eval[2])
        if isfinite(ub_k) && isfinite(ub_v)
            best_ub = min(best_ub, ub_k + ub_v)
        end
        gap = isfinite(best_ub) && !iszero(best_ub) ? (best_ub - current_lb) / abs(best_ub) * 100 : Inf

        oracle.param.dcglp_param.verbose &&
            print_simplex_norm_old_iteration_info(iteration, current_lb, best_ub, gap, ub_k, ub_v, master_time, oracle_times)

        if current_points_feasible
            if current_lb > simplex_norm_old_master_t_value(t_value) + oracle.param.zero_tol
                oracle.param.dcglp_param.verbose &&
                    print_simplex_norm_old_stop_reason("both polar blocks certified in L; generating a disjunctive cut")
                cut = generate_simplex_norm_old_disjunctive_cut(
                    dcglp, current_lb, x_value, length(t_value), zero_indices, one_indices;
                    strengthen = oracle.param.strengthened, lift = oracle.param.lift, zero_tol = oracle.param.zero_tol,
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
                print_simplex_norm_old_stop_reason("both polar blocks are feasible in L; falling back to typical cuts")
            return BendersX.generate_cuts(oracle.typical_oracles[1], x_value, t_value; time_limit = BendersX.get_sec_remaining(start_time, time_limit))
        end

        if gap <= oracle.param.dcglp_param.gap_tolerance
            oracle.param.dcglp_param.verbose &&
                print_simplex_norm_old_stop_reason("dcglp gap tolerance reached")
            break
        end

        if no_improvement >= oracle.param.dcglp_param.halt_limit
            oracle.param.dcglp_param.verbose &&
                print_simplex_norm_old_stop_reason("halt limit reached due to insufficient LB improvement")
            break
        end

        if iteration >= oracle.param.dcglp_param.iter_limit
            oracle.param.dcglp_param.verbose &&
                print_simplex_norm_old_stop_reason("iteration limit reached")
            break
        end

        elapsed_time = time() - start_time
        if elapsed_time >= time_limit || elapsed_time >= oracle.param.dcglp_param.time_limit
            oracle.param.dcglp_param.verbose &&
                print_simplex_norm_old_stop_reason("time limit reached")
            break
        end

        BendersX.add_constraints(dcglp, :con_benders, violated_cuts)
    end

    if current_lb > simplex_norm_old_master_t_value(t_value) + oracle.param.zero_tol
        oracle.param.dcglp_param.verbose &&
            print_simplex_norm_old_stop_reason("LB sufficient; generating a disjunctive cut")
        cut = generate_simplex_norm_old_disjunctive_cut(
            dcglp, current_lb, x_value, length(t_value), zero_indices, one_indices;
            strengthen = oracle.param.strengthened, lift = oracle.param.lift, zero_tol = oracle.param.zero_tol,
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
        print_simplex_norm_old_stop_reason("LB insufficient; falling back to typical cuts")
    return fallback_typical_or_throw(
        oracle,
        x_value,
        t_value,
        start_time,
        time_limit,
        "SimplexNormOldDCGLP terminated without certifying the polar model.";
        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
    )
end
