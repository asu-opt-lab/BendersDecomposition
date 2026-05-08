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
