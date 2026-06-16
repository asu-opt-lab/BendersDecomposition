"""
    AbstractDisjunctiveNormalization

Normalization object that captures the per-variant behavior of a [`SplitOracle`](@ref).
All split-oracle normalization variants in this module are stored on `SplitOracleParam`.
The normalization owns the variant-specific
configuration (`norm`, `core_point_*`, …) and the dispatch surface listed below.

# Normalization interface
- `build_dcglp` — construct the DCGLP for a normalization. Normalizations with
  distinct linking constraints provide their own method.
- `add_normalization_constraint!` — connect the shared `tau`, `sx`, and `st`
  variables according to the normalization.
- `update_dcglp_for_candidate!` — per-call setup for the current
  candidate point (RHS, direction).
- `initialize_dcglp_state` — normalization-aware state initialization.
- `update_dcglp_reference_t!` — normalizations may override to adjust the
  reference epigraph used by the UB recomputation.
- `record_dcglp_oracle_result!` — normalizations may record subproblem
  results into the state.
- `update_dcglp_upper_bound_and_gap!` — normalization-specific UB/gap rule.
- `build_dcglp_disjunctive_cut` — normalization-specific cut extraction.
- `print_dcglp_iteration_info` — normalization-specific iteration log.
"""
abstract type AbstractDisjunctiveNormalization end

"""
    LpDistanceNormalization(norm, adjust_t_to_fx)

Distance-norm DCGLP normalization for `SplitOracle`.
"""
mutable struct LpDistanceNormalization <: AbstractDisjunctiveNormalization
    norm::AbstractNorm
    adjust_t_to_fx::Bool
end

function LpDistanceNormalization(norm::AbstractNorm = LpNorm(Inf); adjust_t_to_fx::Bool = false)
    return LpDistanceNormalization(norm, adjust_t_to_fx)
end

"""
    ReversePolarNormalization()
    ReversePolarNormalization(core_point_x, core_point_t)

Reverse-polar DCGLP normalization for `SplitOracle`. With no core point, this
uses the vertical reverse-polar direction. With a core point, it uses the
direction from the core point to the current candidate.
"""
mutable struct ReversePolarNormalization <: AbstractDisjunctiveNormalization
    core_point_x::Union{Nothing, Vector{Float64}}
    core_point_t::Union{Nothing, Vector{Float64}}
    last_direction_x::Vector{Float64}
    last_direction_t::Vector{Float64}

    function ReversePolarNormalization()
        return new(nothing, nothing, Float64[], Float64[])
    end

    function ReversePolarNormalization(core_point_x::Vector{Float64}, core_point_t::Vector{Float64})
        isempty(core_point_x) && throw(ArgumentError("`core_point_x` must be non-empty."))
        isempty(core_point_t) && throw(ArgumentError("`core_point_t` must be non-empty."))
        return new(copy(core_point_x), copy(core_point_t), Float64[], Float64[])
    end
end

# Default interface implementations (overridden per normalization where needed).

initialize_dcglp_state(::AbstractDisjunctiveNormalization) = DcglpState()

function update_dcglp_for_candidate!(::AbstractDisjunctiveNormalization, oracle::AbstractSplitOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    set_normalized_rhs.(oracle.dcglp[:cont], t_value)
end

update_dcglp_reference_t!(::AbstractDisjunctiveNormalization, ::AbstractSplitOracle, ::Vector{Float64}, t_value::Vector{Float64}, ::Float64, ::Float64) = copy(t_value)

record_dcglp_oracle_result!(::AbstractDisjunctiveNormalization, ::DcglpState, ::Int, ::Vector{Float64}) = nothing

function normalize_add_benders_cuts_to_master(add_benders_cuts_to_master::Union{Bool, Int})
    return add_benders_cuts_to_master === true ? 1 :
           add_benders_cuts_to_master === false ? 0 :
           add_benders_cuts_to_master in (0, 1, 2) ? add_benders_cuts_to_master :
           throw(ArgumentError("`add_benders_cuts_to_master` must be true, false, or an integer in {0, 1, 2}."))
end

function validate_fraction_of_benders_cuts_to_master(fraction::Float64)
    0.0 < fraction <= 1.0 ||
        throw(ArgumentError("`fraction_of_benders_cuts_to_master` must lie in (0, 1]."))
    return fraction
end

function initialize_disjunctive_cut_storage(master::AbstractMaster)
    return [Vector{Hyperplane}() for _ in 1:master.dim_x],
           Hyperplane[],
           Vector{Tuple{SparseVector{Float64, Int}, Float64}}()
end

"""
    SplitOracleParam

Single parameter container for every split-oracle DCGLP normalization variant.
The variant-specific configuration (`norm`, `core_point_*`, …) lives on the
normalization object; this struct holds the configuration that is shared across
all variants plus a reference to the DCGLP loop settings.
"""
mutable struct SplitOracleParam <: AbstractOracleParam
    dcglp_param::DcglpParam
    normalization::AbstractDisjunctiveNormalization
    split_index_selection_rule::SplitIndexSelectionRule
    disjunctive_cut_append_rule::DisjunctiveCutsAppendRule
    add_benders_cuts_to_master::Int
    fraction_of_benders_cuts_to_master::Float64
    reuse_dcglp::Bool
    strengthened::Bool
    lift::Bool
    zero_tol::Float64

    function SplitOracleParam(
        dcglp_param::DcglpParam,
        normalization::AbstractDisjunctiveNormalization,
        split_index_selection_rule::SplitIndexSelectionRule,
        disjunctive_cut_append_rule::DisjunctiveCutsAppendRule,
        add_benders_cuts_to_master::Union{Bool, Int},
        fraction_of_benders_cuts_to_master::Float64,
        reuse_dcglp::Bool,
        strengthened::Bool,
        lift::Bool,
        zero_tol::Float64,
    )
        return new(
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
        SplitOracleParam(disjunctive_norm_param; dcglp_param = DcglpParam(), kwargs...)

    Construct split-oracle parameters from a normalization parameter object, such as
    `LpDistanceNormalization(LpNorm(Inf))` or `ReversePolarNormalization()`.
    """
    function SplitOracleParam(
        normalization::AbstractDisjunctiveNormalization;
        dcglp_param::DcglpParam = DcglpParam(),
        split_index_selection_rule::SplitIndexSelectionRule = RandomFractional(),
        disjunctive_cut_append_rule::DisjunctiveCutsAppendRule = AllDisjunctiveCuts(),
        add_benders_cuts_to_master::Union{Bool, Int} = 1,
        fraction_of_benders_cuts_to_master::Float64 = 1.0,
        reuse_dcglp::Bool = true,
        strengthened::Bool = true,
        lift::Bool = false,
        zero_tol::Float64 = 1.0e-9,
    )
        return new(
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
end

"""
    build_dcglp(master, param)

Build the DCGLP model for a split-oracle normalization. The fallback builder
uses a shared slack layout with `tau`, `sx`, and `st`, then delegates the
normalization-specific connection to `add_normalization_constraint!`.
"""
function build_dcglp(master::AbstractMaster, param::SplitOracleParam)
    return build_dcglp(master, param, param.normalization)
end

function build_dcglp(
    master::AbstractMaster,
    param::SplitOracleParam,
    normalization::AbstractDisjunctiveNormalization,
)
    dcglp = Model(param.dcglp_param.optimizer)
    @variable(dcglp, omega_0[1:2] >= 0)

    @variable(dcglp, omega_x[1:2, 1:master.dim_x])
    @variable(dcglp, omega_t[1:2, 1:master.dim_t])

    @constraint(dcglp, [i in 1:2], omega_t[i, :] .>= -1.0e6 .* omega_0[i])
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
    SplitOracle

Single struct that backs all split-oracle DCGLP normalization variants. The variant identity
is carried by the `normalization` object stored on `param`.
"""
mutable struct SplitOracle <: AbstractSplitOracle
    param::SplitOracleParam
    dcglp::Model
    typical_oracles::Vector{AbstractTypicalOracle}
    disjunctive_cuts_by_index::Vector{Vector{Hyperplane}}
    disjunctive_cuts::Vector{Hyperplane}
    splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}}

    function SplitOracle(
        param::SplitOracleParam,
        dcglp::Model,
        typical_oracles::Vector{AbstractTypicalOracle},
        disjunctive_cuts_by_index::Vector{Vector{Hyperplane}},
        disjunctive_cuts::Vector{Hyperplane},
        splits::Vector{Tuple{SparseVector{Float64, Int}, Float64}},
    )
        return new(param, dcglp, typical_oracles, disjunctive_cuts_by_index, disjunctive_cuts, splits)
    end

    function SplitOracle(
        master::AbstractMaster,
        typical_oracles::NTuple{2, <:AbstractTypicalOracle};
        param::SplitOracleParam = SplitOracleParam(LpDistanceNormalization()),
    )
        dcglp = build_dcglp(master, param)
        cuts_by_index, cuts, splits = initialize_disjunctive_cut_storage(master)

        return new(
            param,
            dcglp,
            AbstractTypicalOracle[typical_oracles...],
            cuts_by_index,
            cuts,
            splits,
        )
    end
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
