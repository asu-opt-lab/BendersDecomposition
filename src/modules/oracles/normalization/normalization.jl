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
    LpDistanceNormalization(norm, adjust_t_to_fx)

Distance-norm DCGLP normalization for `SplitOracle`.
"""
mutable struct LpDistanceNormalization <: AbstractDisjunctiveNormalization
    norm::AbstractNorm
    adjust_t_to_fx::Bool

    function LpDistanceNormalization(norm::AbstractNorm, adjust_t_to_fx::Bool)
        return new(norm, adjust_t_to_fx)
    end

    function LpDistanceNormalization(norm::AbstractNorm = LpNorm(Inf); adjust_t_to_fx::Bool = false)
        return new(norm, adjust_t_to_fx)
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

update_dcglp_reference_t!(::AbstractDisjunctiveNormalization, ::AbstractSplitOracle, ::Vector{Float64}, t_value::Vector{Float64}, ::Float64, ::Float64) = copy(t_value)
