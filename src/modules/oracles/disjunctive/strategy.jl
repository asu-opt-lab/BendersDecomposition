"""
    AbstractDcglpStrategy

Strategy object that captures the per-variant behavior of a [`DcglpOracle`](@ref).
All four DCGLP variants in this module are instances of `DcglpOracle{S}` for some
`S <: AbstractDcglpStrategy`. The strategy owns the variant-specific
configuration (`norm`, `core_point_*`, …) and the dispatch surface listed below.

# Strategy interface (one method per variant)
- [`build_strategy_dcglp`](@ref) — build the JuMP model from `master` / `param`.
- [`update_dcglp_for_candidate!`](@ref) — per-call setup for the current
  candidate point (RHS, direction).
- [`validate_strategy_specific!`](@ref) — extra construction-time validation.
- [`strategy_label`](@ref) — human-readable name used in error messages.
- [`dcglp_tau_value`](@ref), [`dcglp_sx_value`](@ref),
  [`dcglp_lower_bound`](@ref) — solution reads after each DCGLP solve.
- [`initialize_dcglp_state`](@ref) — strategy-aware state initialization.
- [`update_dcglp_reference_t!`](@ref) — strategies may override to adjust the
  reference epigraph used by the UB recomputation.
- [`record_dcglp_oracle_result!`](@ref) — strategies may record subproblem
  results into the state.
- [`update_dcglp_upper_bound_and_gap!`](@ref) — strategy-specific UB/gap rule.
- [`has_dcglp_disjunctive_cut`](@ref) — strategy-specific cut threshold.
- [`build_dcglp_disjunctive_cut`](@ref) — strategy-specific cut extraction.
- [`print_dcglp_iteration_info`](@ref) — strategy-specific iteration log.
"""
abstract type AbstractDcglpStrategy end

"""
    DistanceNormStrategy(norm, adjust_t_to_fx)

Distance-norm DCGLP strategy. Backs [`DistanceNormOracle`](@ref).
"""
mutable struct DistanceNormStrategy <: AbstractDcglpStrategy
    norm::AbstractNorm
    adjust_t_to_fx::Bool
end

"""
    SimplexNormStrategy()

Simplex-style DCGLP strategy. Backs [`SimplexNormOracle`](@ref).
"""
mutable struct SimplexNormStrategy <: AbstractDcglpStrategy end

"""
    VerticalReversePolarStrategy()

Vertical reverse-polar DCGLP strategy. Backs [`VerticalReversePolarOracle`](@ref).
"""
mutable struct VerticalReversePolarStrategy <: AbstractDcglpStrategy end

"""
    DirectionalPolarStrategy(core_point_x, core_point_t)

Directional reverse-polar DCGLP strategy. Backs [`DirectionalPolarOracle`](@ref).
Stores the core point used to compute the directional cut and the direction
vectors from the most recent `generate_cuts` call (`last_direction_*`).
"""
mutable struct DirectionalPolarStrategy <: AbstractDcglpStrategy
    core_point_x::Vector{Float64}
    core_point_t::Vector{Float64}
    last_direction_x::Vector{Float64}
    last_direction_t::Vector{Float64}

    function DirectionalPolarStrategy(core_point_x::Vector{Float64}, core_point_t::Vector{Float64})
        isempty(core_point_x) && throw(ArgumentError("`core_point_x` must be non-empty."))
        isempty(core_point_t) && throw(ArgumentError("`core_point_t` must be non-empty."))
        new(copy(core_point_x), copy(core_point_t), Float64[], Float64[])
    end
end

# Default interface implementations (overridden per strategy where needed).

strategy_label(s::AbstractDcglpStrategy) = string(nameof(typeof(s)))

validate_strategy_specific!(::AbstractDcglpStrategy, ::AbstractMaster) = nothing

initialize_dcglp_state(::AbstractDcglpStrategy) = DcglpState()

update_dcglp_reference_t!(::AbstractDcglpStrategy, ::AbstractDcglpOracle, ::Vector{Float64}, t_value::Vector{Float64}, ::Float64, ::Float64) = copy(t_value)

dcglp_sx_value(::AbstractDcglpStrategy, ::Model) = Float64[]

dcglp_lower_bound(::AbstractDcglpStrategy, dcglp::Model) = objective_value(dcglp)

record_dcglp_oracle_result!(::AbstractDcglpStrategy, ::DcglpState, ::Int, ::Vector{Float64}) = nothing

fallback_for_unexpected_dcglp_status(::AbstractDcglpStrategy) = true
