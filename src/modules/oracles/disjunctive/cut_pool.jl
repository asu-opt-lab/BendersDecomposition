function add_disjunctive_cuts!(oracle::AbstractDcglpOracle, rule::DisjunctiveCutsAppendRule)
    throw(UndefError("update add_disjunctive_cuts! for $(typeof(rule))"))
end

function add_disjunctive_cuts!(oracle::AbstractDcglpOracle, ::NoDisjunctiveCuts)
    delete_registered_constraints!(oracle.dcglp, :con_disjunctive)
end

function add_disjunctive_cuts!(oracle::AbstractDcglpOracle, ::AllDisjunctiveCuts)
    install_disjunctive_cuts!(oracle, oracle.disjunctive_cuts)
end

function add_disjunctive_cuts!(oracle::AbstractDcglpOracle, ::DisjunctiveCutsSmallerIndices)
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

function install_disjunctive_cuts!(oracle::AbstractDcglpOracle, cuts::Vector{Hyperplane})
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

function update_dynamic_dcglp_constraints!(oracle::AbstractDcglpOracle)
    if !oracle.param.reuse_dcglp
        delete_registered_constraints!(oracle.dcglp, :con_benders)
        delete_registered_constraints!(oracle.dcglp, :con_disjunctive)
    end
    add_disjunctive_cuts!(oracle, oracle.param.disjunctive_cut_append_rule)
end

function append_current_disjunctive_cut!(oracle::AbstractDcglpOracle, cut::Hyperplane)
    push!(oracle.disjunctive_cuts, cut)
    if oracle.param.split_index_selection_rule isa SimpleSplit
        index = get_split_index(oracle)
        push!(oracle.disjunctive_cuts_by_index[index], cut)
    end
end

function store_dcglp_disjunctive_cut!(
    oracle::AbstractDcglpOracle,
    cut::Hyperplane,
    master_hyperplanes::Vector{Hyperplane},
    include_disjunctive_cuts_to_hyperplanes::Bool,
)
    append_current_disjunctive_cut!(oracle, cut)
    include_disjunctive_cuts_to_hyperplanes && push!(master_hyperplanes, cut)
end
