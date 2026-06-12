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
    build_dcglp(normalization, master, param)

Build the DCGLP model shared by every split-oracle normalization. Variant
behavior is supplied only through `add_normalization_constraint!`.
"""
function build_dcglp(
    normalization::S,
    master::AbstractMaster,
    param::SplitOracleParam{S},
) where {S <: AbstractDisjunctiveNormalization}
    dcglp = Model(param.dcglp_param.optimizer)
    build_dcglp_skeleton!(dcglp, master)
    add_normalization_constraint!(dcglp, normalization, master, param)
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

function SplitOracle{S}(
    master::AbstractMaster,
    typical_oracles::Vector{T},
    param::SplitOracleParam{S},
) where {S <: AbstractDisjunctiveNormalization, T <: AbstractTypicalOracle}
    label = normalization_label(param.normalization)
    validate_two_typical_oracles!(typical_oracles, label)
    validate_binary_master!(master, label)
    validate_normalization_specific!(param.normalization, master)

    dcglp = build_dcglp(param.normalization, master, param)
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

function SplitOracle(
    master::AbstractMaster,
    typical_oracles::Vector{T},
    param::SplitOracleParam{S},
) where {S <: AbstractDisjunctiveNormalization, T <: AbstractTypicalOracle}
    return SplitOracle{S}(master, typical_oracles, param)
end
