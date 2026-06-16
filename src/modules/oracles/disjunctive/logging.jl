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

function print_dcglp_iteration_info(::LpDistanceNormalization, state::DcglpState, log::DcglpLog)
    print_iteration_info(state, log)
end

function print_dcglp_iteration_info(normalization::ReversePolarNormalization, state::DcglpState, log::DcglpLog)
    if has_core_point(normalization)
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
