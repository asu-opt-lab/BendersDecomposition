"""
    SplitOracleParam{S<:AbstractDisjunctiveNormalizationStrategy}

Single parameter container for every split-oracle DCGLP normalization variant. The variant-specific
configuration (`norm`, `core_point_*`, …) lives on the strategy object; this
struct holds the configuration that is shared across all variants plus a
reference to the DCGLP loop settings.

For read ergonomics the strategy fields are transparently forwarded via
`getproperty`, so existing code can still read `param.norm` or
`param.core_point_x`. Writes are not forwarded: per the package policy
split oracle parameters are fixed at construction. Mutating strategy state
that must stay in sync with the DCGLP (e.g. the directional core point)
goes through dedicated APIs such as [`set_core_point!`](@ref).
"""
mutable struct SplitOracleParam{S<:AbstractDisjunctiveNormalizationStrategy} <: AbstractOracleParam
    dcglp_param::DcglpParam
    strategy::S
    split_index_selection_rule::SplitIndexSelectionRule
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule
    add_benders_cuts_to_master::Int
    fraction_of_benders_cuts_to_master::Float64
    reuse_dcglp::Bool
    strengthened::Bool
    lift::Bool
    zero_tol::Float64
end

"""Parameter container for [`DistanceNormOracle`](@ref)."""
const DistanceNormOracleParam         = SplitOracleParam{DistanceNormStrategy}
"""Parameter container for [`SimplexNormOracle`](@ref)."""
const SimplexNormOracleParam          = SplitOracleParam{SimplexNormStrategy}
"""Parameter container for [`VerticalReversePolarOracle`](@ref)."""
const VerticalReversePolarOracleParam = SplitOracleParam{VerticalReversePolarStrategy}
"""Parameter container for [`DirectionalPolarOracle`](@ref)."""
const DirectionalPolarOracleParam     = SplitOracleParam{DirectionalPolarStrategy}

function Base.getproperty(p::SplitOracleParam, name::Symbol)
    name in fieldnames(SplitOracleParam) ?
        getfield(p, name) :
        getproperty(getfield(p, :strategy), name)
end

function Base.propertynames(p::SplitOracleParam, private::Bool = false)
    return (fieldnames(SplitOracleParam)..., propertynames(getfield(p, :strategy), private)...)
end

"""
    DistanceNormOracleParam(dcglp_param; kwargs...)

Construct the parameters for a [`DistanceNormOracle`](@ref).
"""
function SplitOracleParam{DistanceNormStrategy}(
    dcglp_param::DcglpParam;
    norm::AbstractNorm = LpNorm(Inf),
    adjust_t_to_fx::Bool = false,
    split_index_selection_rule::SplitIndexSelectionRule = RandomFractional(),
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule = AllDisjunctiveCuts(),
    strengthened::Bool = true,
    add_benders_cuts_to_master::Union{Bool, Int} = 1,
    fraction_of_benders_cuts_to_master::Float64 = 1.0,
    reuse_dcglp::Bool = true,
    lift::Bool = false,
    zero_tol::Float64 = 1.0e-9,
)
    return SplitOracleParam(
        dcglp_param,
        DistanceNormStrategy(norm, adjust_t_to_fx),
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
    SimplexNormOracleParam(dcglp_param; kwargs...)

Construct the parameters for a [`SimplexNormOracle`](@ref).
"""
function SplitOracleParam{SimplexNormStrategy}(
    dcglp_param::DcglpParam;
    split_index_selection_rule::SplitIndexSelectionRule = RandomFractional(),
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule = AllDisjunctiveCuts(),
    add_benders_cuts_to_master::Union{Bool, Int} = 1,
    fraction_of_benders_cuts_to_master::Float64 = 1.0,
    reuse_dcglp::Bool = true,
    strengthened::Bool = true,
    lift::Bool = false,
    zero_tol::Float64 = 1.0e-9,
)
    return SplitOracleParam(
        dcglp_param,
        SimplexNormStrategy(),
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
    VerticalReversePolarOracleParam(dcglp_param; kwargs...)

Construct the parameters for a [`VerticalReversePolarOracle`](@ref).
"""
function SplitOracleParam{VerticalReversePolarStrategy}(
    dcglp_param::DcglpParam;
    split_index_selection_rule::SplitIndexSelectionRule = RandomFractional(),
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule = AllDisjunctiveCuts(),
    add_benders_cuts_to_master::Union{Bool, Int} = 1,
    fraction_of_benders_cuts_to_master::Float64 = 1.0,
    reuse_dcglp::Bool = true,
    strengthened::Bool = true,
    lift::Bool = false,
    zero_tol::Float64 = 1.0e-9,
)
    return SplitOracleParam(
        dcglp_param,
        VerticalReversePolarStrategy(),
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
    DirectionalPolarOracleParam(dcglp_param, core_point_x, core_point_t; kwargs...)

Construct the parameters for a [`DirectionalPolarOracle`](@ref).
"""
function SplitOracleParam{DirectionalPolarStrategy}(
    dcglp_param::DcglpParam,
    core_point_x::Vector{Float64},
    core_point_t::Vector{Float64};
    split_index_selection_rule::SplitIndexSelectionRule = RandomFractional(),
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule = AllDisjunctiveCuts(),
    add_benders_cuts_to_master::Union{Bool, Int} = 1,
    fraction_of_benders_cuts_to_master::Float64 = 1.0,
    reuse_dcglp::Bool = true,
    strengthened::Bool = true,
    lift::Bool = false,
    zero_tol::Float64 = 1.0e-9,
)
    return SplitOracleParam(
        dcglp_param,
        DirectionalPolarStrategy(core_point_x, core_point_t),
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
    SplitOracle{S<:AbstractDisjunctiveNormalizationStrategy}

Single struct that backs all four split-oracle DCGLP normalization variants. The variant identity
is carried by the `S` type parameter; the four user-facing names
[`DistanceNormOracle`](@ref), [`SimplexNormOracle`](@ref),
[`VerticalReversePolarOracle`](@ref), [`DirectionalPolarOracle`](@ref) are
`const` aliases for the corresponding parameterization.
"""
mutable struct SplitOracle{S<:AbstractDisjunctiveNormalizationStrategy} <: AbstractSplitOracle
    param::SplitOracleParam{S}
    dcglp::Model
    typical_oracles::Vector{AbstractTypicalOracle}
    disjunctive_cuts_by_index::Vector{Vector{Hyperplane}}
    disjunctive_cuts::Vector{Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}
end

"""Split oracle with distance-norm normalization (`tau ≥ ‖·‖_p`)."""
const DistanceNormOracle        = SplitOracle{DistanceNormStrategy}
"""Split oracle with simplex-style objective `sum(tau)`."""
const SimplexNormOracle         = SplitOracle{SimplexNormStrategy}
"""Vertical reverse-polar split oracle (scalar `tau`, componentwise UB)."""
const VerticalReversePolarOracle = SplitOracle{VerticalReversePolarStrategy}
"""Directional reverse-polar split oracle with configurable core point."""
const DirectionalPolarOracle    = SplitOracle{DirectionalPolarStrategy}

function SplitOracle{S}(
    master::AbstractMaster,
    typical_oracles::Vector{T},
    param::SplitOracleParam{S},
) where {S <: AbstractDisjunctiveNormalizationStrategy, T <: AbstractTypicalOracle}
    label = strategy_label(param.strategy)
    validate_two_typical_oracles!(typical_oracles, label)
    validate_binary_master!(master, label)
    validate_strategy_specific!(param.strategy, master)

    dcglp = build_strategy_dcglp(param.strategy, master, param)
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
) where {S <: AbstractDisjunctiveNormalizationStrategy, T <: AbstractTypicalOracle}
    return SplitOracle{S}(master, typical_oracles, param)
end
