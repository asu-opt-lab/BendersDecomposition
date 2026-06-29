# Utilities shared by SplitOracle.
# Load order follows the dependencies between helper methods.

# -----------------------------------------------------------------------------
# Shared DCGLP helpers
# -----------------------------------------------------------------------------

function delete_registered_constraints!(model::Model, sym::Symbol)
    haskey(model, sym) || return nothing
    registered = model[sym]
    if registered isa AbstractArray
        delete.(Ref(model), registered)
    else
        delete(model, registered)
    end
    unregister(model, sym)
end

function fallback_to_typical_after_dcglp!(
    oracle::AbstractDisjunctiveOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
)
    oracle isa AbstractSplitOracle && rollback_current_dcglp_candidate!(oracle)
    return generate_cuts(
        oracle.typical_oracles[1],
        x_value,
        t_value;
        time_limit = get_sec_remaining(start_time, time_limit),
    )
end

function fallback_typical_or_throw(
    oracle::AbstractDisjunctiveOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
    msg::String;
    throw_typical_cuts_for_errors::Bool,
)
    if throw_typical_cuts_for_errors
        @warn msg
        return fallback_to_typical_after_dcglp!(oracle, x_value, t_value, start_time, time_limit)
    end
    oracle isa AbstractSplitOracle && rollback_current_dcglp_candidate!(oracle)
    throw(UnexpectedModelStatusException(msg))
end

# -----------------------------------------------------------------------------
# Split selection and split constraints
# -----------------------------------------------------------------------------

function select_disjunctive_inequality(x_value::Vector{Float64}, split_selection_rule::SplitIndexSelectionRule; zero_tol = 1.0e-2)
    throw(UndefError("update select_disjunctive_inequality for $(typeof(split_selection_rule))"))
end

function select_disjunctive_inequality(x_value::Vector{Float64}, ::LargestFractional; zero_tol = 1.0e-9)
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(eachindex(x_value))) : maximum(frac_indices)

    phi = spzeros(length(x_value))
    phi[index] = 1.0
    return phi, 0.0
end

function select_disjunctive_inequality(x_value::Vector{Float64}, ::MostFractional; zero_tol = 1.0e-9)
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(eachindex(x_value))) : frac_indices[argmin(abs.(x_value[frac_indices] .- 0.5))]

    phi = spzeros(length(x_value))
    phi[index] = 1.0
    return phi, 0.0
end

function select_disjunctive_inequality(x_value::Vector{Float64}, ::RandomFractional; zero_tol = 1.0e-9)
    frac_indices = filter(i -> zero_tol <= x_value[i] <= 1.0 - zero_tol, eachindex(x_value))
    index = isempty(frac_indices) ? rand(collect(eachindex(x_value))) : rand(frac_indices)

    phi = spzeros(length(x_value))
    phi[index] = 1.0
    return phi, 0.0
end

function get_split_index(oracle::AbstractSplitOracle)
    oracle.param.split_index_selection_rule isa SimpleSplit ||
        throw(AlgorithmException("get_split_index is only valid for simple split rules."))
    isempty(oracle.splits) &&
        throw(AlgorithmException("get_split_index requires at least one selected split."))
    return findfirst(x -> x > 0.5, oracle.splits[end][1])
end

function replace_disjunctive_inequality!(oracle::AbstractSplitOracle)
    dcglp = oracle.dcglp
    phi, phi_0 = oracle.splits[end]

    delete_registered_constraints!(dcglp, :con_split_kappa)
    delete_registered_constraints!(dcglp, :con_split_nu)

    @constraint(dcglp, con_split_kappa, 0 >= dcglp[:omega_0][1] * (phi_0 + 1.0) - phi' * dcglp[:omega_x][1, :])
    @constraint(dcglp, con_split_nu, 0 >= -dcglp[:omega_0][2] * phi_0 + phi' * dcglp[:omega_x][2, :])
end

function retrieve_zero_one(x_value::Vector{Float64}, zero_tol)
    zero_indices = findall(x -> isapprox(x, 0.0; atol = zero_tol), x_value)
    one_indices = findall(x -> isapprox(x, 1.0; atol = zero_tol), x_value)
    return zero_indices, one_indices
end

function add_lifting_constraints!(dcglp::Model, zero_indices::Vector{Int}, one_indices::Vector{Int})
    delete_registered_constraints!(dcglp, :con_zeta)
    delete_registered_constraints!(dcglp, :con_xi)

    !isempty(zero_indices) &&
        @constraint(dcglp, con_zeta[i in 1:2, j in eachindex(zero_indices)], 0 >= dcglp[:omega_x][i, zero_indices[j]])
    !isempty(one_indices) &&
        @constraint(dcglp, con_xi[i in 1:2, j in eachindex(one_indices)], 0 >= dcglp[:omega_0][i] - dcglp[:omega_x][i, one_indices[j]])
end

function choose_split_and_update_lifting!(
    oracle::AbstractSplitOracle,
    x_value::Vector{Float64},
)
    push!(
        oracle.splits,
        select_disjunctive_inequality(
            x_value,
            oracle.param.split_index_selection_rule;
            zero_tol = oracle.param.zero_tol,
        ),
    )
    replace_disjunctive_inequality!(oracle)

    zero_indices, one_indices = oracle.param.lift ? retrieve_zero_one(x_value, oracle.param.zero_tol) : (Int[], Int[])
    add_lifting_constraints!(oracle.dcglp, zero_indices, one_indices)
    return zero_indices, one_indices
end

function rollback_current_dcglp_candidate!(oracle::AbstractSplitOracle)
    !isempty(oracle.splits) && pop!(oracle.splits)
    delete_registered_constraints!(oracle.dcglp, :con_split_kappa)
    delete_registered_constraints!(oracle.dcglp, :con_split_nu)
    delete_registered_constraints!(oracle.dcglp, :con_zeta)
    delete_registered_constraints!(oracle.dcglp, :con_xi)
    delete_registered_constraints!(oracle.dcglp, :initial_L)
    return nothing
end

# -----------------------------------------------------------------------------
# Disjunctive cut pool management
# -----------------------------------------------------------------------------

function add_disjunctive_cuts!(oracle::AbstractSplitOracle, rule::DisjunctiveCutsAppendRule)
    throw(UndefError("update add_disjunctive_cuts! for $(typeof(rule))"))
end

function add_disjunctive_cuts!(oracle::AbstractSplitOracle, ::NoDisjunctiveCuts)
    delete_registered_constraints!(oracle.dcglp, :con_disjunctive)
end

function add_disjunctive_cuts!(oracle::AbstractSplitOracle, ::AllDisjunctiveCuts)
    install_disjunctive_cuts!(oracle, oracle.disjunctive_cuts)
end

function add_disjunctive_cuts!(oracle::AbstractSplitOracle, ::DisjunctiveCutsSmallerIndices)
    oracle.param.split_index_selection_rule isa SimpleSplit ||
        throw(AlgorithmException("DisjunctiveCutsSmallerIndices requires a simple split rule."))

    index = get_split_index(oracle)
    cuts = if index > 1
        reduce(vcat, oracle.disjunctive_cuts_by_index[1:index-1]; init = Hyperplane[])
    else
        Hyperplane[]
    end
    install_disjunctive_cuts!(oracle, cuts)
end

function install_disjunctive_cuts!(oracle::AbstractSplitOracle, cuts::Vector{Hyperplane})
    dcglp = oracle.dcglp
    delete_registered_constraints!(dcglp, :con_disjunctive)
    isempty(cuts) && return nothing

    exprs = AffExpr[]
    for k in 1:2
        append!(
            exprs,
            hyperplanes_to_expression(
                dcglp,
                cuts,
                dcglp[:omega_x][k, :],
                dcglp[:omega_t][k, :],
                dcglp[:omega_0][k],
            ),
        )
    end
    add_constraints(dcglp, :con_disjunctive, exprs)
end

function update_dynamic_dcglp_constraints!(oracle::AbstractSplitOracle)
    if !oracle.param.reuse_dcglp
        delete_registered_constraints!(oracle.dcglp, :con_benders)
        delete_registered_constraints!(oracle.dcglp, :con_disjunctive)
    end
    add_disjunctive_cuts!(oracle, oracle.param.disjunctive_cut_append_rule)
end

function append_current_disjunctive_cut!(oracle::AbstractSplitOracle, cut::Hyperplane)
    push!(oracle.disjunctive_cuts, cut)
    if oracle.param.split_index_selection_rule isa SimpleSplit
        index = get_split_index(oracle)
        push!(oracle.disjunctive_cuts_by_index[index], cut)
    end
end

function store_dcglp_disjunctive_cut!(
    oracle::AbstractSplitOracle,
    cut::Hyperplane,
    master_hyperplanes::Vector{Hyperplane},
    include_disjunctive_cuts_to_hyperplanes::Bool,
)
    append_current_disjunctive_cut!(oracle, cut)
    include_disjunctive_cuts_to_hyperplanes && push!(master_hyperplanes, cut)
end

# -----------------------------------------------------------------------------
# Cut lifting, strengthening, and normalization
# -----------------------------------------------------------------------------

"""
    apply_lift_or_strengthen(dcglp, gamma_x, zero_indices, one_indices; lift, strengthen, zero_tol, gamma_0=0.0)

Apply lifting (preferred if `lift` and any zero/one index exists) or pure
strengthening to the cut coefficients. Returns the (possibly modified)
`(gamma_x, gamma_0)` pair.
"""
function apply_lift_or_strengthen(
    dcglp::Model,
    gamma_x::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int};
    lift::Bool,
    strengthen::Bool,
    zero_tol::Float64,
    gamma_0::Float64 = 0.0,
)
    if lift && (!isempty(zero_indices) || !isempty(one_indices))
        return lift_x_coefficients(
            dcglp, gamma_x, gamma_0, zero_indices, one_indices;
            strengthen = strengthen, zero_tol = zero_tol,
        )
    elseif strengthen
        sigma, delta = read_strengthening_duals(dcglp)
        return -strengthen_coefficients(-gamma_x, sigma, delta; zero_tol = zero_tol), gamma_0
    end
    return gamma_x, gamma_0
end

read_strengthening_duals(dcglp::Model) = (
    Dict(1 => dual(dcglp[:con_split_kappa]), 2 => dual(dcglp[:con_split_nu])),
    Dict(1 => dual.(dcglp[:condelta][1, :]), 2 => dual.(dcglp[:condelta][2, :])),
)

"""
    lift_x_coefficients(dcglp, gamma_x, gamma_0, zero_indices, one_indices; strengthen, zero_tol)

Apply lifting (and optional strengthening) to `gamma_x`. Returns the lifted
`(gamma_x, gamma_0)` pair in the common cut orientation.
"""
function lift_x_coefficients(
    dcglp::Model,
    gamma_x::Vector{Float64},
    gamma_0::Float64,
    zero_indices::Vector{Int},
    one_indices::Vector{Int};
    strengthen::Bool,
    zero_tol::Float64,
)
    zeta_k = !isempty(zero_indices) ? dual.(dcglp[:con_zeta][1, :]) : Float64[]
    zeta_v = !isempty(zero_indices) ? dual.(dcglp[:con_zeta][2, :]) : Float64[]
    xi_k = !isempty(one_indices) ? dual.(dcglp[:con_xi][1, :]) : Float64[]
    xi_v = !isempty(one_indices) ? dual.(dcglp[:con_xi][2, :]) : Float64[]

    lifted_gamma_0 = gamma_0 - sum(max.(xi_k, xi_v))
    lifted_gamma_x = -gamma_x
    lifted_gamma_x[zero_indices] .= -gamma_x[zero_indices] .+ max.(zeta_k, zeta_v)
    lifted_gamma_x[one_indices] .= -gamma_x[one_indices] .- max.(xi_k, xi_v)

    if strengthen
        sigma, base_delta = read_strengthening_duals(dcglp)
        delta_1 = copy(base_delta[1])
        delta_2 = copy(base_delta[2])
        delta_1[zero_indices] .+= -zeta_k .+ max.(zeta_k, zeta_v)
        delta_2[zero_indices] .+= -zeta_v .+ max.(zeta_k, zeta_v)
        lifted_gamma_x = strengthen_coefficients(lifted_gamma_x, sigma, Dict(1 => delta_1, 2 => delta_2); zero_tol = zero_tol)
    end

    return -lifted_gamma_x, lifted_gamma_0
end

function strengthen_coefficients(gamma_x, sigma, delta; zero_tol = 1.0e-9)
    a_1 = gamma_x .- delta[1]
    a_2 = gamma_x .- delta[2]
    sigma_sum = sigma[1] + sigma[2]
    sigma_sum < zero_tol && return copy(gamma_x)
    m = (a_1 .- a_2) ./ sigma_sum
    return min.(a_1 .- sigma[1] .* floor.(m), a_2 .+ sigma[2] .* ceil.(m))
end

# -----------------------------------------------------------------------------
# DCGLP logging
# -----------------------------------------------------------------------------

function format_sparse_terms(coeffs::SparseVector{Float64, Int}, var_name::String; zero_tol::Float64 = 1.0e-10)
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

function format_hyperplane(h::Hyperplane; zero_tol::Float64 = 1.0e-10)
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

function print_dcglp_iteration_info(::AbstractDisjunctiveNormalization, state::DcglpState, log::DcglpLog)
    print_iteration_info(state, log)
end

function print_disjunctive_cut(
    oracle::SplitOracle,
    cut::Hyperplane,
    x_value::Vector{Float64},
    t_value::Vector{Float64};
    zero_tol::Float64 = 1.0e-10,
)
    println("   Disjunctive cut:")
    println("      " * format_hyperplane(cut; zero_tol = zero_tol))
    if oracle.param.split_index_selection_rule isa SimpleSplit
        @printf("   Split index: x[%d]\n", get_split_index(oracle))
    end
    @printf("   Cut violation at current point: %.6f\n", evaluate_violation(cut, x_value, t_value))
end

# -----------------------------------------------------------------------------
# DCGLP solve loop
# -----------------------------------------------------------------------------

"""
    should_fallback_typical(normalization, oracle, x_value, t_value) -> Bool

If `true`, `generate_cuts` short-circuits and delegates to the first typical
oracle without solving the DCGLP. Default implementation returns `false`.
The directional normalization uses this hook to skip the DCGLP when the
candidate point coincides with the core point.
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
    if should_fallback_typical(normalization, oracle, x_value, reference_t)
        return fallback_to_typical_after_dcglp!(oracle, x_value, t_value, start_time, time_limit)
    end

    while true
        state = DcglpState()
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
                    rollback_current_dcglp_candidate!(oracle)
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
        cut = build_dcglp_disjunctive_cut(
            normalization,
            dcglp,
            oracle.param,
            zero_indices,
            one_indices,
        )
        oracle.param.dcglp_param.verbose && print_disjunctive_cut(oracle, cut, x_value, t_value; zero_tol = oracle.param.zero_tol)
        store_dcglp_disjunctive_cut!(oracle, cut, hyperplanes, include_disjunctive_cuts_to_hyperplanes)
        return false, hyperplanes, fill(Inf, length(t_value))
    end

    return fallback_to_typical_after_dcglp!(oracle, x_value, t_value, start_time, time_limit)
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

function disjunctive_cut_normalization_value(
    normalization::AbstractDisjunctiveNormalization,
    dcglp::Model,
    gamma_x::Vector{Float64},
    gamma_t::Vector{Float64},
)
    throw(UndefError("update disjunctive_cut_normalization_value for $(typeof(normalization))"))
end

function build_dcglp_disjunctive_cut(
    normalization::AbstractDisjunctiveNormalization,
    dcglp::Model,
    common::SplitOracleParam,
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

    norm_value = disjunctive_cut_normalization_value(
        normalization,
        dcglp,
        gamma_x,
        gamma_t,
    )
    return Hyperplane(gamma_x ./ norm_value, gamma_t ./ norm_value, gamma_0 / norm_value)
end
