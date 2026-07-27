# -----------------------------------------------------------------------------
# Reverse-polar normalization
# -----------------------------------------------------------------------------

"""
    ReversePolarNormalization()
    ReversePolarNormalization(core_point_x, core_point_t)
    ReversePolarNormalization(; core_point_x, core_point_t)
    ReversePolarNormalization(; core_direction_x, core_direction_t)

Reverse-polar DCGLP normalization for `SplitOracle`. With no core point or
direction, this uses the vertical reverse-polar direction. With a core point,
it uses the direction from the core point to the current candidate. With a core
direction, it uses that fixed direction for every candidate.
"""
_reverse_polar_vector(::Nothing, ::String) = nothing
_reverse_polar_vector(vector::AbstractVector{<:Real}, ::String) = Float64.(collect(vector))

mutable struct ReversePolarNormalization <: AbstractDisjunctiveNormalization
    core_point_x::Union{Nothing, Vector{Float64}}
    core_point_t::Union{Nothing, Vector{Float64}}
    core_direction_x::Union{Nothing, Vector{Float64}}
    core_direction_t::Union{Nothing, Vector{Float64}}
    use_core_point::Bool

    function ReversePolarNormalization(
        core_point_x::Union{Nothing, AbstractVector{<:Real}},
        core_point_t::Union{Nothing, AbstractVector{<:Real}},
        core_direction_x::Union{Nothing, AbstractVector{<:Real}},
        core_direction_t::Union{Nothing, AbstractVector{<:Real}},
    )
        core_point_x = _reverse_polar_vector(core_point_x, "core_point_x")
        core_point_t = _reverse_polar_vector(core_point_t, "core_point_t")
        core_direction_x = _reverse_polar_vector(core_direction_x, "core_direction_x")
        core_direction_t = _reverse_polar_vector(core_direction_t, "core_direction_t")

        has_point = core_point_x !== nothing || core_point_t !== nothing
        has_direction = core_direction_x !== nothing || core_direction_t !== nothing
        use_core_point = has_point

        has_point && has_direction &&
            throw(ArgumentError("Provide either `core_point_x`/`core_point_t` or `core_direction_x`/`core_direction_t`, not both."))

        if has_point
            core_point_x !== nothing ||
                throw(ArgumentError("`core_point_x` must be provided when `core_point_t` is provided."))
            core_point_t !== nothing ||
                throw(ArgumentError("`core_point_t` must be provided when `core_point_x` is provided."))
            isempty(core_point_x) && throw(ArgumentError("`core_point_x` must be non-empty."))
            isempty(core_point_t) && throw(ArgumentError("`core_point_t` must be non-empty."))
        end

        if has_direction
            core_direction_x !== nothing ||
                throw(ArgumentError("`core_direction_x` must be provided when `core_direction_t` is provided."))
            core_direction_t !== nothing ||
                throw(ArgumentError("`core_direction_t` must be provided when `core_direction_x` is provided."))
            isempty(core_direction_x) && throw(ArgumentError("`core_direction_x` must be non-empty."))
            isempty(core_direction_t) && throw(ArgumentError("`core_direction_t` must be non-empty."))
            LinearAlgebra.norm(vcat(core_direction_x, core_direction_t), Inf) > 0.0 ||
                throw(ArgumentError("`core_direction_x`/`core_direction_t` must not define the zero direction."))
        end

        return new(
            core_point_x,
            core_point_t,
            core_direction_x,
            core_direction_t,
            use_core_point,
        )
    end

    function ReversePolarNormalization(core_point_x::AbstractVector{<:Real}, core_point_t::AbstractVector{<:Real})
        return ReversePolarNormalization(core_point_x, core_point_t, nothing, nothing)
    end

    function ReversePolarNormalization(;
        core_point_x::Union{Nothing, AbstractVector{<:Real}} = nothing,
        core_point_t::Union{Nothing, AbstractVector{<:Real}} = nothing,
        core_direction_x::Union{Nothing, AbstractVector{<:Real}} = nothing,
        core_direction_t::Union{Nothing, AbstractVector{<:Real}} = nothing,
    )
        return ReversePolarNormalization(core_point_x, core_point_t, core_direction_x, core_direction_t)
    end
end

function check_reverse_polar_dimension(vector::Vector{Float64}, expected::Int, name::String)
    length(vector) == expected ||
        throw(DimensionMismatch("`$name` has length $(length(vector)) but expected $expected."))
    return nothing
end

function prepare_disjunctive_normalization!(normalization::ReversePolarNormalization, master::AbstractMaster)
    if normalization.use_core_point
        check_reverse_polar_dimension(normalization.core_point_x, master.dim_x, "core_point_x")
        check_reverse_polar_dimension(normalization.core_point_t, master.dim_t, "core_point_t")
        return nothing
    end

    if normalization.core_direction_x === nothing && normalization.core_direction_t === nothing
        normalization.core_direction_x = zeros(master.dim_x)
        normalization.core_direction_t = ones(master.dim_t)
        return nothing
    end

    check_reverse_polar_dimension(normalization.core_direction_x, master.dim_x, "core_direction_x")
    check_reverse_polar_dimension(normalization.core_direction_t, master.dim_t, "core_direction_t")
    return nothing
end

function add_normalization_constraint!(
    normalization::ReversePolarNormalization,
    dcglp::Model,
    tau::VariableRef,
    sx::AbstractVector{VariableRef},
    st::AbstractVector{VariableRef},
)
    @constraint(dcglp, con_reverse_polar_x[j in eachindex(sx)], sx[j] + 0.0 * tau == 0.0)
    @constraint(dcglp, con_reverse_polar_t[j in eachindex(st)], st[j] + 0.0 * tau == 0.0)
    return nothing
end

function should_fallback_typical(normalization::ReversePolarNormalization, oracle::SplitOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    normalization.use_core_point || return false

    direction_x = x_value .- normalization.core_point_x
    direction_t = t_value .- normalization.core_point_t

    if LinearAlgebra.norm(vcat(direction_x, direction_t), Inf) <= oracle.param.zero_tol
        return true
    end

    return false
end

function update_dcglp_for_candidate!(normalization::ReversePolarNormalization, oracle::SplitOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    set_normalized_rhs.(oracle.dcglp[:cont], t_value)

    direction_x, direction_t =
        normalization.use_core_point ?
        (x_value .- normalization.core_point_x, t_value .- normalization.core_point_t) :
        (.-normalization.core_direction_x, .-normalization.core_direction_t)
    update_reverse_polar_constraints!(oracle, direction_x, direction_t)
    return nothing
end

function update_dcglp_reference_t!(
    normalization::ReversePolarNormalization,
    oracle::SplitOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    start_time::Float64,
    time_limit::Float64,
)
    oracle.param.adjust_t_to_fx || return copy(t_value)

    reference_t = adjust_dcglp_reference_t_to_fx!(oracle, x_value, t_value, start_time, time_limit)
    direction_x, direction_t =
        normalization.use_core_point ?
        (x_value .- normalization.core_point_x, reference_t .- normalization.core_point_t) :
        (.-normalization.core_direction_x, .-normalization.core_direction_t)
    update_reverse_polar_constraints!(oracle, direction_x, direction_t)
    return reference_t
end

function update_dcglp_upper_bound_and_gap!(
    normalization::ReversePolarNormalization,
    state::DcglpState,
    log::DcglpLog,
    reference_t::Vector{Float64},
    t_value::Vector{Float64},
)
    # todo: consider computing upper bound based on projection onto half-ray; Currently it's based on `is_in_L` indicator.

    fill_dcglp_omega_t_estimates!(state, t_value)
    previous_ub = log.n_iter >= 1 ? log.iterations[end].UB : Inf
    state.UB = state.is_in_L[1] && state.is_in_L[2] ? min(previous_ub, state.LB) : previous_ub
    state.gap =
        isfinite(state.UB) ?
        max(state.UB - state.LB, 0.0) / max(abs(state.UB), 1.0) * 100 :
        Inf
    return nothing
end

function disjunctive_cut_normalization_value(
    normalization::ReversePolarNormalization,
    dcglp::Model,
    gamma_x::Vector{Float64},
    gamma_t::Vector{Float64},
)
    tau = dcglp[:tau]
    direction_x = [normalized_coefficient(con, tau) for con in dcglp[:con_reverse_polar_x]]
    direction_t = [normalized_coefficient(con, tau) for con in dcglp[:con_reverse_polar_t]]
    direction_value = dot(gamma_x, direction_x) + dot(gamma_t, direction_t)
    !iszero(direction_value) ||
        throw(AlgorithmException("ReversePolarNormalization cut normalization failed because the directional support is numerically zero."))
    return direction_value
end

"""
    set_core_point!(oracle::SplitOracle, core_point_x, core_point_t)

Update the core point used by a directional reverse-polar `SplitOracle`. The supplied
vectors must have the same dimensions as the oracle's existing core point.
"""
function set_core_point!(oracle::SplitOracle, core_point_x::Vector{Float64}, core_point_t::Vector{Float64})
    return set_core_point!(oracle.param.normalization, core_point_x, core_point_t)
end

function set_core_point!(normalization::ReversePolarNormalization, core_point_x::Vector{Float64}, core_point_t::Vector{Float64})
    normalization.use_core_point ||
        throw(ArgumentError("`set_core_point!` requires a core-point ReversePolarNormalization."))
    length(core_point_x) == length(normalization.core_point_x) ||
        throw(DimensionMismatch("`core_point_x` has length $(length(core_point_x)) but expected $(length(normalization.core_point_x))."))
    length(core_point_t) == length(normalization.core_point_t) ||
        throw(DimensionMismatch("`core_point_t` has length $(length(core_point_t)) but expected $(length(normalization.core_point_t))."))

    normalization.core_point_x .= core_point_x
    normalization.core_point_t .= core_point_t
    return nothing
end

function update_reverse_polar_constraints!(
    oracle::SplitOracle,
    direction_x::Vector{Float64},
    direction_t::Vector{Float64},
)
    dcglp = oracle.dcglp
    tau = dcglp[:tau]

    for j in eachindex(direction_x)
        set_normalized_coefficient(dcglp[:con_reverse_polar_x][j], tau, direction_x[j])
    end
    for j in eachindex(direction_t)
        set_normalized_coefficient(dcglp[:con_reverse_polar_t][j], tau, direction_t[j])
    end

    return nothing
end
