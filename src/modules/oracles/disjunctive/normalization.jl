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
- `normalization_label` — human-readable name used in error messages.
- `dcglp_sx_value` — normalization-specific solution reads after each DCGLP solve.
- `initialize_dcglp_state` — normalization-aware state initialization.
- `update_dcglp_reference_t!` — normalizations may override to adjust the
  reference epigraph used by the UB recomputation.
- `record_dcglp_oracle_result!` — normalizations may record subproblem
  results into the state.
- `update_dcglp_upper_bound_and_gap!` — normalization-specific UB/gap rule.
- `has_dcglp_disjunctive_cut` — cut threshold rule.
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

normalization_label(s::AbstractDisjunctiveNormalization) = string(nameof(typeof(s)))

initialize_dcglp_state(::AbstractDisjunctiveNormalization) = DcglpState()

update_dcglp_reference_t!(::AbstractDisjunctiveNormalization, ::AbstractSplitOracle, ::Vector{Float64}, t_value::Vector{Float64}, ::Float64, ::Float64) = copy(t_value)

dcglp_sx_value(::AbstractDisjunctiveNormalization, ::Model) = Float64[]

record_dcglp_oracle_result!(::AbstractDisjunctiveNormalization, ::DcglpState, ::Int, ::Vector{Float64}) = nothing

has_dcglp_disjunctive_cut(::AbstractDisjunctiveNormalization, current_lb::Float64, ::Vector{Float64}, zero_tol::Float64) = current_lb >= zero_tol

fallback_for_unexpected_dcglp_status(::AbstractDisjunctiveNormalization) = true
