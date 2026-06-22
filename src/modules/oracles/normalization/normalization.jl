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
- `update_dcglp_reference_t!` — normalizations may override to adjust the
  reference epigraph used by the UB recomputation.
- `update_dcglp_upper_bound_and_gap!` — normalization-specific UB/gap rule.
- `build_dcglp_disjunctive_cut` — normalization-specific cut extraction.
- `print_dcglp_iteration_info` — normalization-specific iteration log.
"""
abstract type AbstractDisjunctiveNormalization end

"""
    LpDistanceNormalization(norm)

Distance-norm DCGLP normalization for `SplitOracle`.
"""
mutable struct LpDistanceNormalization <: AbstractDisjunctiveNormalization
    norm_p::AbstractNorm

    function LpDistanceNormalization(norm::AbstractNorm = LpNorm(Inf))
        return new(norm)
    end
end

"""
    ReversePolarNormalization()
    ReversePolarNormalization(core_point_x, core_point_t)
    ReversePolarNormalization(; core_point_x, core_point_t)
    ReversePolarNormalization(; core_diretion_x, core_diretion_t)

Reverse-polar DCGLP normalization for `SplitOracle`. With no core point or
direction, this uses the vertical reverse-polar direction. With a core point,
it uses the direction from the core point to the current candidate. With a core
direction, it uses that fixed direction for every candidate.
"""
mutable struct ReversePolarNormalization <: AbstractDisjunctiveNormalization
    core_point_x::Union{Nothing, Vector{Float64}}
    core_point_t::Union{Nothing, Vector{Float64}}
    core_diretion_x::Union{Nothing, Vector{Float64}}
    core_diretion_t::Union{Nothing, Vector{Float64}}

    function ReversePolarNormalization(
        core_point_x::Union{Nothing, Vector{Float64}},
        core_point_t::Union{Nothing, Vector{Float64}},
        core_diretion_x::Union{Nothing, Vector{Float64}},
        core_diretion_t::Union{Nothing, Vector{Float64}},
    )
        has_point = core_point_x !== nothing || core_point_t !== nothing
        has_direction = core_diretion_x !== nothing || core_diretion_t !== nothing
        default_vertical_direction = !has_point && !has_direction

        has_point && has_direction &&
            throw(ArgumentError("Provide either `core_point_x`/`core_point_t` or `core_diretion_x`/`core_diretion_t`, not both."))

        if default_vertical_direction
            core_diretion_x = Float64[]
            core_diretion_t = [1.0]
            has_direction = true
        end

        if has_point
            core_point_x !== nothing ||
                throw(ArgumentError("`core_point_x` must be provided when `core_point_t` is provided."))
            core_point_t !== nothing ||
                throw(ArgumentError("`core_point_t` must be provided when `core_point_x` is provided."))
            isempty(core_point_x) && throw(ArgumentError("`core_point_x` must be non-empty."))
            isempty(core_point_t) && throw(ArgumentError("`core_point_t` must be non-empty."))
        end

        if has_direction
            core_diretion_x !== nothing ||
                throw(ArgumentError("`core_diretion_x` must be provided when `core_diretion_t` is provided."))
            core_diretion_t !== nothing ||
                throw(ArgumentError("`core_diretion_t` must be provided when `core_diretion_x` is provided."))
            !default_vertical_direction && isempty(core_diretion_x) &&
                throw(ArgumentError("`core_diretion_x` must be non-empty."))
            isempty(core_diretion_t) && throw(ArgumentError("`core_diretion_t` must be non-empty."))
            LinearAlgebra.norm(vcat(core_diretion_x, core_diretion_t), Inf) > 0.0 ||
                throw(ArgumentError("`core_diretion_x`/`core_diretion_t` must not define the zero direction."))
        end

        return new(
            core_point_x === nothing ? nothing : copy(core_point_x),
            core_point_t === nothing ? nothing : copy(core_point_t),
            core_diretion_x === nothing ? nothing : copy(core_diretion_x),
            core_diretion_t === nothing ? nothing : copy(core_diretion_t),
        )
    end

    function ReversePolarNormalization(core_point_x::Vector{Float64}, core_point_t::Vector{Float64})
        return ReversePolarNormalization(core_point_x, core_point_t, nothing, nothing)
    end

    function ReversePolarNormalization(;
        core_point_x::Union{Nothing, Vector{Float64}} = nothing,
        core_point_t::Union{Nothing, Vector{Float64}} = nothing,
        core_diretion_x::Union{Nothing, Vector{Float64}} = nothing,
        core_diretion_t::Union{Nothing, Vector{Float64}} = nothing,
    )
        return ReversePolarNormalization(core_point_x, core_point_t, core_diretion_x, core_diretion_t)
    end
end

# Default interface implementations (overridden per normalization where needed).

function update_dcglp_for_candidate!(::AbstractDisjunctiveNormalization, oracle::AbstractSplitOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    set_normalized_rhs.(oracle.dcglp[:cont], t_value)
end

function adjust_dcglp_reference_t_to_fx!(
    oracle::AbstractSplitOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
)
    dcglp = oracle.dcglp
    delete_registered_constraints!(dcglp, :initial_L)

    _, initial_hyperplanes, f_x = generate_cuts(
        oracle.typical_oracles[1],
        x_value,
        t_value;
        time_limit = get_sec_remaining(start_time, time_limit),
    )
    any(isnan, f_x) &&
        throw(AlgorithmException("solve_dcglp!: `t_value` cannot be adjusted to `f(x)` since $(typeof(oracle.typical_oracles[1])) does not compute `f(x)`."))

    initial_benders_cuts = AffExpr[]
    for k in 1:2
        append!(
            initial_benders_cuts,
            hyperplanes_to_expression(
                dcglp,
                initial_hyperplanes,
                dcglp[:omega_x][k, :],
                dcglp[:omega_t][k, :],
                dcglp[:omega_0][k],
            ),
        )
    end
    dcglp[:initial_L] = @constraint(dcglp, 0 .>= initial_benders_cuts)
    set_normalized_rhs.(dcglp[:cont], f_x)
    return copy(f_x)
end

function update_dcglp_reference_t!(
    ::AbstractDisjunctiveNormalization,
    oracle::AbstractSplitOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
)
    oracle.param.adjust_t_to_fx || return copy(t_value)
    return adjust_dcglp_reference_t_to_fx!(oracle, x_value, t_value, start_time, time_limit)
end
