function install_disjunctive_cuts!(oracle::SimplexNormOldDCGLPOracle, cuts::Vector{BendersX.Hyperplane})
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

function append_current_disjunctive_cut!(oracle::SimplexNormOldDCGLPOracle, cut::BendersX.Hyperplane)
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

function add_disjunctive_cuts!(oracle::SimplexNormOldDCGLPOracle, ::BendersX.NoDisjunctiveCuts)
    return nothing
end

function add_disjunctive_cuts!(oracle::SimplexNormOldDCGLPOracle, ::BendersX.AllDisjunctiveCuts)
    return nothing
end

function add_disjunctive_cuts!(oracle::SimplexNormOldDCGLPOracle, ::BendersX.DisjunctiveCutsSmallerIndices)
    current_index = get_split_index(oracle)
    reusable = BendersX.Hyperplane[]
    for idx in 1:(current_index - 1)
        append!(reusable, oracle.disjunctiveCutsByIndex[idx])
    end
    install_disjunctive_cuts!(oracle, reusable)
    return nothing
end

function BendersX.generate_cuts(
    oracle::SimplexNormOldDCGLPOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64};
    time_limit::Float64 = 3600.0,
    throw_typical_cuts_for_errors::Bool = true,
    include_disjunctive_cuts_to_hyperplanes::Bool = true,
)
    push!(
        oracle.splits,
        BendersX.select_disjunctive_inequality(x_value, oracle.param.split_index_selection_rule; zero_tol = oracle.param.zero_tol),
    )
    replace_disjunctive_inequality!(oracle)

    if !oracle.param.reuse_dcglp
        delete_registered_constraints!(oracle.dcglp, :con_benders)
    end
    add_disjunctive_cuts!(oracle, oracle.param.disjunctive_cut_append_rule)

    JuMP.set_normalized_rhs.(oracle.dcglp[:conx], x_value)

    zero_indices, one_indices = oracle.param.lift ? BendersX.retrieve_zero_one(x_value, oracle.param.zero_tol) : (Int[], Int[])
    BendersX.add_lifting_constraints!(oracle.dcglp, zero_indices, one_indices)

    start_time = time()

    return optimize_simplex_norm_old_dcglp!(
        oracle,
        x_value,
        t_value,
        zero_indices,
        one_indices;
        start_time = start_time,
        time_limit = time_limit,
        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
        include_disjunctive_cuts_to_hyperplanes = include_disjunctive_cuts_to_hyperplanes,
    )
end
