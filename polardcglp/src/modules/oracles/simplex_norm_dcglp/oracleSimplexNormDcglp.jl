mutable struct SimplexNormDCGLPParam <: BendersX.AbstractOracleParam
    dcglp_param::BendersX.DcglpParam
    split_index_selection_rule::BendersX.SplitIndexSelectionRule
    disjunctive_cut_append_rule::BendersX.DisjunctiveCutsAppendRule
    add_benders_cuts_to_master::Int
    fraction_of_benders_cuts_to_master::Float64
    reuse_dcglp::Bool
    strengthened::Bool
    lift::Bool
    zero_tol::Float64

    function SimplexNormDCGLPParam(
        dcglp_param::BendersX.DcglpParam;
        split_index_selection_rule::BendersX.SplitIndexSelectionRule = BendersX.RandomFractional(),
        disjunctive_cut_append_rule::BendersX.DisjunctiveCutsAppendRule = BendersX.AllDisjunctiveCuts(),
        add_benders_cuts_to_master::Union{Bool, Int} = 1,
        fraction_of_benders_cuts_to_master::Float64 = 1.0,
        reuse_dcglp::Bool = true,
        strengthened::Bool = true,
        lift::Bool = false,
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
            strengthened,
            lift,
            zero_tol,
        )
    end
end

mutable struct SimplexNormDCGLPOracle <: BendersX.AbstractDisjunctiveOracle
    param::SimplexNormDCGLPParam
    dcglp::Model
    typical_oracles::Vector{BendersX.AbstractTypicalOracle}
    disjunctiveCutsByIndex::Vector{Vector{BendersX.Hyperplane}}
    disjunctiveCuts::Vector{BendersX.Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}

    function SimplexNormDCGLPOracle(
        master::BendersX.AbstractMaster,
        typical_oracles::Vector{T},
        param::SimplexNormDCGLPParam,
    ) where {T <: BendersX.AbstractTypicalOracle}
        length(typical_oracles) == 2 ||
            throw(ArgumentError("SimplexNormDCGLPOracle requires exactly two typical oracles."))

        for xi in master.x
            is_binary(xi) || throw(ArgumentError("SimplexNormDCGLPOracle requires all master variables to be binary."))
        end

        dcglp = build_simplex_norm_dcglp(master, param)
        disjunctive_cuts_by_index = [Vector{BendersX.Hyperplane}() for _ in 1:master.dim_x]
        splits = Vector{Tuple{SparseVector{Float64, Int}, Float64}}()

        new(param, dcglp, Vector{BendersX.AbstractTypicalOracle}(typical_oracles), disjunctive_cuts_by_index, BendersX.Hyperplane[], splits)
    end
end

function BendersX.generate_cuts(
    oracle::SimplexNormDCGLPOracle,
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

    return BendersX.solve_dcglp!(
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

function build_simplex_norm_dcglp(master::BendersX.AbstractMaster, param::SimplexNormDCGLPParam)
    dcglp = Model(param.dcglp_param.optimizer)

    @variable(dcglp, tau[1:master.dim_t])
    @variable(dcglp, omega_0[1:2] >= 0)
    @variable(dcglp, omega_x[1:2, 1:master.dim_x])
    @variable(dcglp, omega_t[1:2, 1:master.dim_t])

    @objective(dcglp, Min, sum(tau))

    @constraint(dcglp, [i in 1:2], omega_t[i, :] .>= SIMPLEXNORM_T_LOWER_BOUND .* omega_0[i])
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

function get_split_index(oracle::SimplexNormDCGLPOracle)
    isa(oracle.param.split_index_selection_rule, BendersX.SimpleSplit) ||
        throw(BendersX.AlgorithmException("get_split_index is only valid for simple split rules."))
    return findfirst(x -> x > 0.5, oracle.splits[end][1])
end

function replace_disjunctive_inequality!(oracle::SimplexNormDCGLPOracle)
    dcglp = oracle.dcglp
    phi, phi_0 = oracle.splits[end]

    delete_registered_constraints!(dcglp, :con_split_kappa)
    delete_registered_constraints!(dcglp, :con_split_nu)

    @constraint(dcglp, con_split_kappa, 0 >= dcglp[:omega_0][1] * (phi_0 + 1.0) - phi' * dcglp[:omega_x][1, :])
    @constraint(dcglp, con_split_nu, 0 >= -dcglp[:omega_0][2] * phi_0 + phi' * dcglp[:omega_x][2, :])

    return nothing
end
