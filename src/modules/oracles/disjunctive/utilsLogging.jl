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
