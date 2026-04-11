module PolarDCGLP

using JuMP
using MathOptInterface
using LinearAlgebra
using Printf
using Random
using SparseArrays

import BendersX

const MOI = MathOptInterface
const POLAR_T_LOWER_BOUND = -1e6

export PolarDCGLPParam, PolarDCGLPOracle

mutable struct PolarDCGLPParam <: BendersX.AbstractOracleParam
    dcglp_param::BendersX.DcglpParam
    split_index_selection_rule::BendersX.SplitIndexSelectionRule
    disjunctive_cut_append_rule::BendersX.DisjunctiveCutsAppendRule
    add_benders_cuts_to_master::Int
    fraction_of_benders_cuts_to_master::Float64
    reuse_dcglp::Bool
    adjust_t_to_fx::Bool
    zero_tol::Float64

    function PolarDCGLPParam(
        dcglp_param::BendersX.DcglpParam;
        split_index_selection_rule::BendersX.SplitIndexSelectionRule = BendersX.RandomFractional(),
        disjunctive_cut_append_rule::BendersX.DisjunctiveCutsAppendRule = BendersX.AllDisjunctiveCuts(),
        add_benders_cuts_to_master::Union{Bool, Int} = 1,
        fraction_of_benders_cuts_to_master::Float64 = 1.0,
        reuse_dcglp::Bool = true,
        adjust_t_to_fx::Bool = false,
        zero_tol::Float64 = 1e-9,
    )
        add_bcuts_to_master =
            add_benders_cuts_to_master === true ? 1 :
            add_benders_cuts_to_master === false ? 0 :
            add_benders_cuts_to_master in (0, 1, 2) ? add_benders_cuts_to_master :
            throw(ArgumentError("`add_benders_cuts_to_master` must be true, false, or an integer in {0, 1, 2}."))

        0.0 < fraction_of_benders_cuts_to_master <= 1.0 ||
            throw(ArgumentError("`fraction_of_benders_cuts_to_master` must lie in (0, 1]."))

        new(
            dcglp_param,
            split_index_selection_rule,
            disjunctive_cut_append_rule,
            add_bcuts_to_master,
            fraction_of_benders_cuts_to_master,
            reuse_dcglp,
            adjust_t_to_fx,
            zero_tol,
        )
    end
end

mutable struct PolarDCGLPOracle <: BendersX.AbstractDisjunctiveOracle
    param::PolarDCGLPParam
    dcglp::Model
    typical_oracles::Vector{BendersX.AbstractTypicalOracle}
    disjunctiveCutsByIndex::Vector{Vector{BendersX.Hyperplane}}
    disjunctiveCuts::Vector{BendersX.Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}

    function PolarDCGLPOracle(
        master::BendersX.AbstractMaster,
        typical_oracles::Vector{T},
        param::PolarDCGLPParam,
    ) where {T <: BendersX.AbstractTypicalOracle}
        length(typical_oracles) == 2 ||
            throw(ArgumentError("PolarDCGLPOracle requires exactly two typical oracles."))

        for xi in master.x
            is_binary(xi) || throw(ArgumentError("PolarDCGLPOracle requires all master variables to be binary."))
        end

        dcglp = build_polar_dcglp(master, param)
        disjunctive_cuts_by_index = [Vector{BendersX.Hyperplane}() for _ in 1:master.dim_x]
        splits = Vector{Tuple{SparseVector{Float64, Int}, Float64}}()

        new(param, dcglp, Vector{BendersX.AbstractTypicalOracle}(typical_oracles), disjunctive_cuts_by_index, BendersX.Hyperplane[], splits)
    end
end

function build_polar_dcglp(master::BendersX.AbstractMaster, param::PolarDCGLPParam)
    dcglp = Model(param.dcglp_param.optimizer)

    @variable(dcglp, tau[1:master.dim_t])
    @variable(dcglp, omega_0[1:2] >= 0)
    @variable(dcglp, omega_x[1:2, 1:master.dim_x])
    @variable(dcglp, omega_t[1:2, 1:master.dim_t])

    @objective(dcglp, Min, sum(tau))

    @constraint(dcglp, [i in 1:2], omega_t[i, :] .>= POLAR_T_LOWER_BOUND .* omega_0[i])
    @constraint(dcglp, coneta[i in 1:2, j in 1:master.dim_x], 0 >= -omega_0[i] + omega_x[i, j])
    @constraint(dcglp, condelta[i in 1:2, j in 1:master.dim_x], 0 >= -omega_x[i, j])

    @constraint(dcglp, con0, omega_0[1] + omega_0[2] == 1.0)
    @constraint(dcglp, conx[j in 1:master.dim_x], omega_x[1, j] + omega_x[2, j] == 0.0)
    @constraint(dcglp, cont[j in 1:master.dim_t], omega_t[1, j] + omega_t[2, j] - tau[j] == 0.0)

    for i in 1:2
        BendersX.transfer_scaled_linear_rows_and_bounds_with_types!(master.model, master.x, dcglp, dcglp[:omega_x][i, :], dcglp[:omega_0][i])
    end

    return dcglp
end

function delete_registered_constraints!(model::Model, sym::Symbol)
    haskey(model, sym) || return
    registered = model[sym]
    if registered isa AbstractArray
        delete.(model, registered)
    else
        delete(model, registered)
    end
    unregister(model, sym)
end

remaining_time(start_time::Float64, time_limit::Float64; tol::Float64 = 1e-4) =
    max(time_limit - (time() - start_time), tol)

function print_polar_iteration_info(
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

function update_polar_no_improvement!(
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

function print_polar_stop_reason(msg::String)
    @printf("   Stop: %s\n", msg)
    return nothing
end

function split_phi_and_rhs(x_value::Vector{Float64}, ::BendersX.LargestFractional; zero_tol::Float64 = 1e-9)
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(eachindex(x_value))) : maximum(frac_indices)

    phi = spzeros(Float64, length(x_value))
    phi[index] = 1.0
    return phi, 0.0
end

function split_phi_and_rhs(x_value::Vector{Float64}, ::BendersX.MostFractional; zero_tol::Float64 = 1e-9)
    gap_x = abs.(x_value .- 0.5)
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(eachindex(x_value))) : argmin(gap_x)

    phi = spzeros(Float64, length(x_value))
    phi[index] = 1.0
    return phi, 0.0
end

function split_phi_and_rhs(x_value::Vector{Float64}, ::BendersX.RandomFractional; zero_tol::Float64 = 1e-9)
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(eachindex(x_value))) : rand(frac_indices)

    phi = spzeros(Float64, length(x_value))
    phi[index] = 1.0
    return phi, 0.0
end

function split_phi_and_rhs(
    x_value::Vector{Float64},
    rule::BendersX.SplitIndexSelectionRule;
    zero_tol::Float64 = 1e-9,
)
    throw(BendersX.UndefError("PolarDCGLP does not implement split selection for $(typeof(rule))."))
end

function get_split_index(oracle::PolarDCGLPOracle)
    isa(oracle.param.split_index_selection_rule, BendersX.SimpleSplit) ||
        throw(BendersX.AlgorithmException("get_split_index is only valid for simple split rules."))
    return findfirst(x -> x > 0.5, oracle.splits[end][1])
end

function replace_disjunctive_inequality!(oracle::PolarDCGLPOracle)
    dcglp = oracle.dcglp
    phi, phi_0 = oracle.splits[end]

    delete_registered_constraints!(dcglp, :con_split_kappa)
    delete_registered_constraints!(dcglp, :con_split_nu)

    @constraint(dcglp, con_split_kappa, 0 >= dcglp[:omega_0][1] * (phi_0 + 1.0) - phi' * dcglp[:omega_x][1, :])
    @constraint(dcglp, con_split_nu, 0 >= -dcglp[:omega_0][2] * phi_0 + phi' * dcglp[:omega_x][2, :])

    return nothing
end

function install_disjunctive_cuts!(oracle::PolarDCGLPOracle, cuts::Vector{BendersX.Hyperplane})
    dcglp = oracle.dcglp
    delete_registered_constraints!(dcglp, :con_disjunctive)
    isempty(cuts) && return nothing

    exprs = AffExpr[]
    for i in 1:2
        append!(exprs, BendersX.hyperplanes_to_expression(dcglp, cuts, dcglp[:omega_x][i, :], dcglp[:omega_t][i, :], dcglp[:omega_0][i]))
    end

    BendersX.add_constraints(dcglp, :con_disjunctive, exprs)
    return nothing
end

function add_disjunctive_cuts!(oracle::PolarDCGLPOracle, ::BendersX.NoDisjunctiveCuts)
    return nothing
end

function add_disjunctive_cuts!(oracle::PolarDCGLPOracle, ::BendersX.AllDisjunctiveCuts)
    return nothing
end

function add_disjunctive_cuts!(oracle::PolarDCGLPOracle, ::BendersX.DisjunctiveCutsSmallerIndices)
    current_index = get_split_index(oracle)
    reusable = BendersX.Hyperplane[]
    for idx in 1:(current_index - 1)
        append!(reusable, oracle.disjunctiveCutsByIndex[idx])
    end
    install_disjunctive_cuts!(oracle, reusable)
    return nothing
end

function append_current_disjunctive_cut!(oracle::PolarDCGLPOracle, cut::BendersX.Hyperplane)
    push!(oracle.disjunctiveCuts, cut)
    if isa(oracle.param.split_index_selection_rule, BendersX.SimpleSplit)
        push!(oracle.disjunctiveCutsByIndex[get_split_index(oracle)], cut)
    end
    if isa(oracle.param.disjunctive_cut_append_rule, BendersX.AllDisjunctiveCuts)
        dcglp = oracle.dcglp
        exprs = AffExpr[]
        for i in 1:2
            append!(exprs, BendersX.hyperplanes_to_expression(dcglp, [cut], dcglp[:omega_x][i, :], dcglp[:omega_t][i, :], dcglp[:omega_0][i]))
        end
        BendersX.add_constraints(dcglp, :con_disjunctive, exprs)
    end
    return nothing
end

function fallback_typical_or_throw(
    oracle::PolarDCGLPOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
    msg::String;
    throw_typical_cuts_for_errors::Bool,
)
    if throw_typical_cuts_for_errors
        @warn msg
        return BendersX.generate_cuts(oracle.typical_oracles[1], x_value, t_value; time_limit = remaining_time(start_time, time_limit))
    end
    throw(BendersX.UnexpectedModelStatusException(msg))
end

hyperplane_violation(h::BendersX.Hyperplane, x_value::Vector{Float64}, t_value::Vector{Float64}) =
    h.a_0 + dot(h.a_x, x_value) + dot(h.a_t, t_value)

polar_master_t_value(t_value::Vector{Float64}) = sum(t_value)

function seed_initial_l!(oracle::PolarDCGLPOracle, x_value::Vector{Float64}, t_value::Vector{Float64}, start_time::Float64, time_limit::Float64)
    dcglp = oracle.dcglp
    delete_registered_constraints!(dcglp, :initial_L)

    t_seed = copy(t_value)
    if oracle.param.adjust_t_to_fx
        _, _, f_x = BendersX.generate_cuts(oracle.typical_oracles[1], x_value, t_value; time_limit = remaining_time(start_time, time_limit))
        any(isnan, f_x) && throw(BendersX.AlgorithmException("PolarDCGLP cannot adjust `t_value` to `f(x)` because the typical oracle returned NaN."))
        t_seed = f_x
    end

    _, seed_cuts, _ = BendersX.generate_cuts(oracle.typical_oracles[1], x_value, t_seed; time_limit = remaining_time(start_time, time_limit))
    exprs = AffExpr[]
    for i in 1:2
        append!(exprs, BendersX.hyperplanes_to_expression(dcglp, seed_cuts, dcglp[:omega_x][i, :], dcglp[:omega_t][i, :], dcglp[:omega_0][i]))
    end
    dcglp[:initial_L] = @constraint(dcglp, 0 .>= exprs)
    return nothing
end

function optimize_polar_dcglp!(
    oracle::PolarDCGLPOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64};
    start_time::Float64,
    time_limit::Float64,
    throw_typical_cuts_for_errors::Bool,
    include_disjunctive_cuts_to_hyperplanes::Bool,
)
    dcglp = oracle.dcglp
    master_hyperplanes = BendersX.Hyperplane[]

    best_lb = -Inf
    best_ub = Inf
    no_improvement = 0
    iteration = 0
    prev_lb = -Inf

    while true
        iteration += 1
        JuMP.set_time_limit_sec(dcglp, remaining_time(start_time, time_limit))

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
                "PolarDCGLP master failed during optimization: $(err)";
                throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
            )
        end

        status = termination_status(dcglp)
        if status == TIME_LIMIT
            throw(BendersX.TimeLimitException("Time limit reached during PolarDCGLP solving."))
        elseif status != OPTIMAL
            return fallback_typical_or_throw(
                oracle,
                x_value,
                t_value,
                start_time,
                time_limit,
                "PolarDCGLP master terminated with unexpected status $(status).";
                throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
            )
        end

        current_lb = objective_value(dcglp)
        if iteration >= 2
            no_improvement = update_polar_no_improvement!(
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
                    time_limit = remaining_time(start_time, time_limit),
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
            print_polar_iteration_info(iteration, current_lb, best_ub, gap, ub_k, ub_v, master_time, oracle_times)

        if current_points_feasible
            if current_lb > polar_master_t_value(t_value) + oracle.param.zero_tol
                oracle.param.dcglp_param.verbose &&
                    print_polar_stop_reason("both polar blocks certified in L; generating a disjunctive cut")
                gamma_x = dual.(dcglp[:conx])
                gamma_0 = current_lb - dot(gamma_x, x_value)
                cut = BendersX.Hyperplane(gamma_x, fill(-1.0, length(t_value)), gamma_0)
                append_current_disjunctive_cut!(oracle, cut)
                if include_disjunctive_cuts_to_hyperplanes
                    push!(master_hyperplanes, cut)
                end
                return false, master_hyperplanes, fill(Inf, length(t_value))
            end

            oracle.param.dcglp_param.verbose &&
                print_polar_stop_reason("both polar blocks are feasible in L; falling back to typical cuts")
            return BendersX.generate_cuts(oracle.typical_oracles[1], x_value, t_value; time_limit = remaining_time(start_time, time_limit))
        end

        if gap <= oracle.param.dcglp_param.gap_tolerance
            oracle.param.dcglp_param.verbose &&
                print_polar_stop_reason("dcglp gap tolerance reached before full certification; falling back to typical cuts")
            return fallback_typical_or_throw(
                oracle,
                x_value,
                t_value,
                start_time,
                time_limit,
                "PolarDCGLP reached the dcglp gap tolerance before certifying the polar model.";
                throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
            )
        end

        BendersX.add_constraints(dcglp, :con_benders, violated_cuts)

        if no_improvement >= oracle.param.dcglp_param.halt_limit
            oracle.param.dcglp_param.verbose &&
                print_polar_stop_reason("halt limit reached due to insufficient LB improvement; falling back to typical cuts")
            return fallback_typical_or_throw(
                oracle,
                x_value,
                t_value,
                start_time,
                time_limit,
                "PolarDCGLP stalled before certifying the polar model.";
                throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
            )
        end

        if iteration >= oracle.param.dcglp_param.iter_limit
            oracle.param.dcglp_param.verbose &&
                print_polar_stop_reason("iteration limit reached before certification; falling back to typical cuts")
            return fallback_typical_or_throw(
                oracle,
                x_value,
                t_value,
                start_time,
                time_limit,
                "PolarDCGLP reached the iteration limit before certifying the polar model.";
                throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
            )
        end

        elapsed_time = time() - start_time
        if elapsed_time >= time_limit || elapsed_time >= oracle.param.dcglp_param.time_limit
            oracle.param.dcglp_param.verbose &&
                print_polar_stop_reason("time limit reached before certification; falling back to typical cuts")
            return fallback_typical_or_throw(
                oracle,
                x_value,
                t_value,
                start_time,
                time_limit,
                "PolarDCGLP reached the time limit before certifying the polar model.";
                throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
            )
        end
    end
end

function BendersX.generate_cuts(
    oracle::PolarDCGLPOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64};
    time_limit::Float64 = 3600.0,
    throw_typical_cuts_for_errors::Bool = true,
    include_disjunctive_cuts_to_hyperplanes::Bool = true,
)
    push!(
        oracle.splits,
        split_phi_and_rhs(x_value, oracle.param.split_index_selection_rule; zero_tol = oracle.param.zero_tol),
    )
    replace_disjunctive_inequality!(oracle)

    if !oracle.param.reuse_dcglp
        delete_registered_constraints!(oracle.dcglp, :con_benders)
    end
    add_disjunctive_cuts!(oracle, oracle.param.disjunctive_cut_append_rule)

    JuMP.set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    start_time = time()
    seed_initial_l!(oracle, x_value, t_value, start_time, time_limit)

    return optimize_polar_dcglp!(
        oracle,
        x_value,
        t_value;
        start_time = start_time,
        time_limit = time_limit,
        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
        include_disjunctive_cuts_to_hyperplanes = include_disjunctive_cuts_to_hyperplanes,
    )
end

function BendersX.root_node_processing!(master::BendersX.Master, preprocessing::BendersX.DisjunctiveRootNodePreprocessing)
    root_param = deepcopy(preprocessing.params)
    undo = JuMP.relax_integrality(master.model)

    start_time = time()
    try
        env_root_typical = preprocessing.seq_type(master, preprocessing.typical_oracle; param = root_param)
        BendersX.solve!(env_root_typical)

        root_param.time_limit -= time() - start_time

        env_root_disjunctive = preprocessing.seq_type(master, preprocessing.disjunctive_oracle; param = root_param)
        BendersX.solve!(env_root_disjunctive)
    finally
        undo()
    end

    return time() - start_time
end

end # module PolarDCGLP
