"""
    DcglpOracleParam{S<:AbstractDcglpStrategy}

Single parameter container for every DCGLP oracle variant. The variant-specific
configuration (`norm`, `core_point_*`, …) lives on the strategy object; this
struct holds the configuration that is shared across all variants plus a
reference to the DCGLP loop settings.

For read ergonomics the strategy fields are transparently forwarded via
`getproperty`, so existing code can still read `param.norm` or
`param.core_point_x`. Writes are not forwarded: per the package policy
DCGLP oracle parameters are fixed at construction. Mutating strategy state
that must stay in sync with the DCGLP (e.g. the directional core point)
goes through dedicated APIs such as [`set_core_point!`](@ref).
"""
mutable struct DcglpOracleParam{S<:AbstractDcglpStrategy} <: AbstractOracleParam
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
const DistanceNormOracleParam         = DcglpOracleParam{DistanceNormStrategy}
"""Parameter container for [`SimplexNormOracle`](@ref)."""
const SimplexNormOracleParam          = DcglpOracleParam{SimplexNormStrategy}
"""Parameter container for [`VerticalReversePolarOracle`](@ref)."""
const VerticalReversePolarOracleParam = DcglpOracleParam{VerticalReversePolarStrategy}
"""Parameter container for [`DirectionalPolarOracle`](@ref)."""
const DirectionalPolarOracleParam     = DcglpOracleParam{DirectionalPolarStrategy}

function Base.getproperty(p::DcglpOracleParam, name::Symbol)
    name in fieldnames(DcglpOracleParam) ?
        getfield(p, name) :
        getproperty(getfield(p, :strategy), name)
end

function Base.propertynames(p::DcglpOracleParam, private::Bool = false)
    return (fieldnames(DcglpOracleParam)..., propertynames(getfield(p, :strategy), private)...)
end

"""
    DistanceNormOracleParam(dcglp_param; kwargs...)

Construct the parameters for a [`DistanceNormOracle`](@ref).
"""
function DcglpOracleParam{DistanceNormStrategy}(
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
    return DcglpOracleParam(
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
function DcglpOracleParam{SimplexNormStrategy}(
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
    return DcglpOracleParam(
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
function DcglpOracleParam{VerticalReversePolarStrategy}(
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
    return DcglpOracleParam(
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
function DcglpOracleParam{DirectionalPolarStrategy}(
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
    return DcglpOracleParam(
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
    DcglpOracle{S<:AbstractDcglpStrategy}

Single struct that backs all four DCGLP oracle variants. The variant identity
is carried by the `S` type parameter; the four user-facing names
[`DistanceNormOracle`](@ref), [`SimplexNormOracle`](@ref),
[`VerticalReversePolarOracle`](@ref), [`DirectionalPolarOracle`](@ref) are
`const` aliases for the corresponding parameterization.
"""
mutable struct DcglpOracle{S<:AbstractDcglpStrategy} <: AbstractDcglpOracle
    param::DcglpOracleParam{S}
    dcglp::Model
    typical_oracles::Vector{AbstractTypicalOracle}
    disjunctive_cuts_by_index::Vector{Vector{Hyperplane}}
    disjunctive_cuts::Vector{Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}
end

"""DCGLP oracle with distance-norm normalization (`tau ≥ ‖·‖_p`)."""
const DistanceNormOracle        = DcglpOracle{DistanceNormStrategy}
"""DCGLP oracle with simplex-style objective `sum(tau)`."""
const SimplexNormOracle         = DcglpOracle{SimplexNormStrategy}
"""Vertical reverse-polar DCGLP oracle (scalar `tau`, componentwise UB)."""
const VerticalReversePolarOracle = DcglpOracle{VerticalReversePolarStrategy}
"""Directional reverse-polar DCGLP oracle with configurable core point."""
const DirectionalPolarOracle    = DcglpOracle{DirectionalPolarStrategy}

function DcglpOracle{S}(
    master::AbstractMaster,
    typical_oracles::Vector{T},
    param::DcglpOracleParam{S},
) where {S <: AbstractDcglpStrategy, T <: AbstractTypicalOracle}
    label = strategy_label(param.strategy)
    validate_two_typical_oracles!(typical_oracles, label)
    validate_binary_master!(master, label)
    validate_strategy_specific!(param.strategy, master)

    dcglp = build_strategy_dcglp(param.strategy, master, param)
    cuts_by_index, cuts, splits = initialize_disjunctive_cut_storage(master)

    return DcglpOracle{S}(
        param,
        dcglp,
        Vector{AbstractTypicalOracle}(typical_oracles),
        cuts_by_index,
        cuts,
        splits,
    )
end
