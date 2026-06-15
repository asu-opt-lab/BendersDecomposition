"""
    SplitOracleParam{S<:AbstractDisjunctiveNormalization}

Single parameter container for every split-oracle DCGLP normalization variant. The variant-specific
configuration (`norm`, `core_point_*`, …) lives on the normalization object; this
struct holds the configuration that is shared across all variants plus a
reference to the DCGLP loop settings.

For read ergonomics the normalization fields are transparently forwarded via
`getproperty`, so existing code can still read `param.norm` or
`param.core_point_x`. Writes are not forwarded: per the package policy
split oracle parameters are fixed at construction. Mutating normalization state
that must stay in sync with the DCGLP (e.g. the directional core point)
goes through dedicated APIs such as [`set_core_point!`](@ref).
"""
mutable struct SplitOracleParam{S<:AbstractDisjunctiveNormalization} <: AbstractOracleParam
    dcglp_param::DcglpParam
    normalization::S
    split_index_selection_rule::SplitIndexSelectionRule
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule
    add_benders_cuts_to_master::Int
    fraction_of_benders_cuts_to_master::Float64
    reuse_dcglp::Bool
    strengthened::Bool
    lift::Bool
    zero_tol::Float64
end

function Base.getproperty(p::SplitOracleParam, name::Symbol)
    name in fieldnames(SplitOracleParam) ?
        getfield(p, name) :
        getproperty(getfield(p, :normalization), name)
end

function Base.propertynames(p::SplitOracleParam, private::Bool = false)
    return (fieldnames(SplitOracleParam)..., propertynames(getfield(p, :normalization), private)...)
end

"""
    SplitOracleParam(disjunctive_norm_param; dcglp_param = DcglpParam(), kwargs...)

Construct split-oracle parameters from a normalization parameter object, such as
`LpDistanceNormalization(LpNorm(Inf))` or `EpigraphSumNormalization()`.
"""
function SplitOracleParam(
    normalization::S;
    dcglp_param::DcglpParam = DcglpParam(),
    split_index_selection_rule::SplitIndexSelectionRule = RandomFractional(),
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule = AllDisjunctiveCuts(),
    add_benders_cuts_to_master::Union{Bool, Int} = 1,
    fraction_of_benders_cuts_to_master::Float64 = 1.0,
    reuse_dcglp::Bool = true,
    strengthened::Bool = true,
    lift::Bool = false,
    zero_tol::Float64 = 1.0e-9,
) where {S <: AbstractDisjunctiveNormalization}
    return SplitOracleParam(
        dcglp_param,
        normalization,
        split_index_selection_rule,
        disjunctive_cut_append_rule,
        normalize_add_benders_cuts_to_master(add_benders_cuts_to_master),
        validate_fraction_of_benders_cuts_to_master(fraction_of_benders_cuts_to_master),
        reuse_dcglp,
        strengthened,
        lift,
        zero_tol,
    )
end

"""
    build_dcglp(master, param)

Build the DCGLP model for a split-oracle normalization. The fallback builder
uses the shared distance-normalization layout and delegates the normalization
cone to `add_normalization_constraint!`; other normalizations may provide a
more specific method.
"""
function build_dcglp(
    master::AbstractMaster,
    param::SplitOracleParam{S},
) where {S <: AbstractDisjunctiveNormalization}
    normalization = param.normalization
    dcglp = Model(param.dcglp_param.optimizer)
    @variable(dcglp, omega_0[1:2] >= 0)

    @variable(dcglp, omega_x[1:2, 1:master.dim_x])
    @variable(dcglp, omega_t[1:2, 1:master.dim_t])

    @constraint(dcglp, [i in 1:2], omega_t[i, :] .>= DCGLP_OMEGA_T_LOWER_BOUND .* omega_0[i])
    @constraint(dcglp, coneta[i in 1:2, j in 1:master.dim_x], 0 >= -omega_0[i] + omega_x[i, j])
    @constraint(dcglp, condelta[i in 1:2, j in 1:master.dim_x], 0 >= -omega_x[i, j])

    @constraint(dcglp, con0, omega_0[1] + omega_0[2] == 1)

    for i in 1:2
        transfer_scaled_linear_rows_and_bounds_with_types!(
            master.model,
            master.x,
            dcglp,
            omega_x[i, :],
            omega_0[i],
        )
    end

    @variable(dcglp, tau)
    @variable(dcglp, sx[1:master.dim_x])
    @variable(dcglp, st[1:master.dim_t])

    @objective(dcglp, Min, tau)

    @constraint(dcglp, conx, dcglp[:omega_x][1, :] + dcglp[:omega_x][2, :] - sx .== 0)
    @constraint(dcglp, cont[j = 1:master.dim_t], dcglp[:omega_t][1, j] + dcglp[:omega_t][2, j] - st[j] == 0)

    add_normalization_constraint!(normalization, dcglp, tau, sx, st)
    return dcglp
end

"""
    SplitOracle{S<:AbstractDisjunctiveNormalization}

Single struct that backs all four split-oracle DCGLP normalization variants. The variant identity
is carried by the `S` type parameter.
"""
mutable struct SplitOracle{S<:AbstractDisjunctiveNormalization} <: AbstractSplitOracle
    param::SplitOracleParam{S}
    dcglp::Model
    typical_oracles::Vector{AbstractTypicalOracle}
    disjunctive_cuts_by_index::Vector{Vector{Hyperplane}}
    disjunctive_cuts::Vector{Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}
end

function SplitOracle(
    master::AbstractMaster,
    typical_oracles::Vector{T},
    normalization::S;
    kwargs...,
) where {S <: AbstractDisjunctiveNormalization, T <: AbstractTypicalOracle}
    param = SplitOracleParam(normalization; kwargs...)
    label = normalization_label(normalization)
    validate_two_typical_oracles!(typical_oracles, label)
    validate_binary_master!(master, label)
    validate_normalization_specific!(normalization, master)

    dcglp = build_dcglp(master, param)
    cuts_by_index, cuts, splits = initialize_disjunctive_cut_storage(master)

    return SplitOracle{S}(
        param,
        dcglp,
        Vector{AbstractTypicalOracle}(typical_oracles),
        cuts_by_index,
        cuts,
        splits,
    )
end
"""
    generate_cuts(oracle::SplitOracle, x_value, t_value; kwargs...)

Single entry point for every DCGLP oracle variant. Per-variant behavior is
dispatched through the normalization attached to `oracle.param.normalization`.
"""
function generate_cuts(
    oracle::SplitOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64};
    time_limit = 3600.0,
    throw_typical_cuts_for_errors::Bool = true,
    include_disjunctive_cuts_to_hyperplanes::Bool = true,
)
    normalization = oracle.param.normalization

    should_fallback_typical(normalization, oracle, x_value, t_value) &&
        return generate_cuts(oracle.typical_oracles[1], x_value, t_value; time_limit = Float64(time_limit))

    start_time = time()
    zero_indices, one_indices = choose_split_and_update_lifting!(oracle, x_value)
    update_dynamic_dcglp_constraints!(oracle)
    update_dcglp_for_candidate!(normalization, oracle, x_value, t_value)

    return solve_dcglp_loop!(
        oracle,
        x_value,
        t_value,
        zero_indices,
        one_indices;
        start_time = start_time,
        time_limit = Float64(time_limit),
        throw_typical_cuts_for_errors = throw_typical_cuts_for_errors,
        include_disjunctive_cuts_to_hyperplanes = include_disjunctive_cuts_to_hyperplanes,
    )
end

