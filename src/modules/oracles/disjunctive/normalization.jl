"""
    AbstractDisjunctiveNormalization

Normalization object that captures the per-variant behavior of a [`SplitOracle`](@ref).
All four split-oracle normalization variants in this module are instances of `SplitOracle{S}` for some
`S <: AbstractDisjunctiveNormalization`. The normalization owns the variant-specific
configuration (`norm`, `core_point_*`, …) and the dispatch surface listed below.

# Normalization interface (one method per variant)
- [`build_normalization_dcglp`](@ref) — build the JuMP model from `master` / `param`.
- [`update_dcglp_for_candidate!`](@ref) — per-call setup for the current
  candidate point (RHS, direction).
- [`validate_normalization_specific!`](@ref) — extra construction-time validation.
- [`normalization_label`](@ref) — human-readable name used in error messages.
- [`dcglp_tau_value`](@ref), [`dcglp_sx_value`](@ref),
  [`dcglp_lower_bound`](@ref) — solution reads after each DCGLP solve.
- [`initialize_dcglp_state`](@ref) — normalization-aware state initialization.
- [`update_dcglp_reference_t!`](@ref) — normalizations may override to adjust the
  reference epigraph used by the UB recomputation.
- [`record_dcglp_oracle_result!`](@ref) — normalizations may record subproblem
  results into the state.
- [`update_dcglp_upper_bound_and_gap!`](@ref) — normalization-specific UB/gap rule.
- [`has_dcglp_disjunctive_cut`](@ref) — normalization-specific cut threshold.
- [`build_dcglp_disjunctive_cut`](@ref) — normalization-specific cut extraction.
- [`print_dcglp_iteration_info`](@ref) — normalization-specific iteration log.
"""
abstract type AbstractDisjunctiveNormalization end

"""
    LpDistanceNormalization(norm, adjust_t_to_fx)

Distance-norm DCGLP normalization. Backs `SplitOracle{LpDistanceNormalization}`.
"""
mutable struct LpDistanceNormalization <: AbstractDisjunctiveNormalization
    norm::AbstractNorm
    adjust_t_to_fx::Bool
end

"""
    EpigraphSumNormalization()

Epigraph-sum DCGLP normalization. Backs `SplitOracle{EpigraphSumNormalization}`.
"""
mutable struct EpigraphSumNormalization <: AbstractDisjunctiveNormalization end

"""
    VerticalReversePolarNormalization()

Vertical reverse-polar DCGLP normalization. Backs `SplitOracle{VerticalReversePolarNormalization}`.
"""
mutable struct VerticalReversePolarNormalization <: AbstractDisjunctiveNormalization end

"""
    DirectionalReversePolarNormalization(core_point_x, core_point_t)

Directional reverse-polar DCGLP normalization. Backs `SplitOracle{DirectionalReversePolarNormalization}`.
Stores the core point used to compute the directional cut and the direction
vectors from the most recent `generate_cuts` call (`last_direction_*`).
"""
mutable struct DirectionalReversePolarNormalization <: AbstractDisjunctiveNormalization
    core_point_x::Vector{Float64}
    core_point_t::Vector{Float64}
    last_direction_x::Vector{Float64}
    last_direction_t::Vector{Float64}

    function DirectionalReversePolarNormalization(core_point_x::Vector{Float64}, core_point_t::Vector{Float64})
        isempty(core_point_x) && throw(ArgumentError("`core_point_x` must be non-empty."))
        isempty(core_point_t) && throw(ArgumentError("`core_point_t` must be non-empty."))
        new(copy(core_point_x), copy(core_point_t), Float64[], Float64[])
    end
end

# Default interface implementations (overridden per normalization where needed).

normalization_label(s::AbstractDisjunctiveNormalization) = string(nameof(typeof(s)))

validate_normalization_specific!(::AbstractDisjunctiveNormalization, ::AbstractMaster) = nothing

initialize_dcglp_state(::AbstractDisjunctiveNormalization) = DcglpState()

update_dcglp_reference_t!(::AbstractDisjunctiveNormalization, ::AbstractSplitOracle, ::Vector{Float64}, t_value::Vector{Float64}, ::Float64, ::Float64) = copy(t_value)

dcglp_sx_value(::AbstractDisjunctiveNormalization, ::Model) = Float64[]

dcglp_lower_bound(::AbstractDisjunctiveNormalization, dcglp::Model) = objective_value(dcglp)

record_dcglp_oracle_result!(::AbstractDisjunctiveNormalization, ::DcglpState, ::Int, ::Vector{Float64}) = nothing

fallback_for_unexpected_dcglp_status(::AbstractDisjunctiveNormalization) = true
