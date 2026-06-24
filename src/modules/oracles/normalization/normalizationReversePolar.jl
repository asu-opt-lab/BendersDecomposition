# -----------------------------------------------------------------------------
# Reverse-polar normalization
# -----------------------------------------------------------------------------

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

has_core_point(normalization::ReversePolarNormalization) = normalization.core_point_x !== nothing
has_core_direction(normalization::ReversePolarNormalization) = normalization.core_diretion_x !== nothing
is_default_vertical_direction(normalization::ReversePolarNormalization) =
    !has_core_point(normalization) &&
    normalization.core_diretion_x !== nothing &&
    normalization.core_diretion_t !== nothing &&
    isempty(normalization.core_diretion_x) &&
    normalization.core_diretion_t == [1.0]

function check_reverse_polar_dimension(vector::Vector{Float64}, expected::Int, name::String)
    length(vector) == expected ||
        throw(DimensionMismatch("`$name` has length $(length(vector)) but expected $expected."))
    return nothing
end

function prepare_disjunctive_normalization!(normalization::ReversePolarNormalization, master::AbstractMaster)
    if has_core_point(normalization)
        check_reverse_polar_dimension(normalization.core_point_x, master.dim_x, "core_point_x")
        check_reverse_polar_dimension(normalization.core_point_t, master.dim_t, "core_point_t")
        return nothing
    end

    if is_default_vertical_direction(normalization)
        normalization.core_diretion_x = zeros(master.dim_x)
        normalization.core_diretion_t = ones(master.dim_t)
        return nothing
    end

    if has_core_direction(normalization)
        check_reverse_polar_dimension(normalization.core_diretion_x, master.dim_x, "core_diretion_x")
        check_reverse_polar_dimension(normalization.core_diretion_t, master.dim_t, "core_diretion_t")
        return nothing
    end

    throw(ArgumentError("ReversePolarNormalization requires a core point, a core direction, or the default vertical direction."))
end

function reverse_polar_direction(normalization::ReversePolarNormalization, x_value::Vector{Float64}, t_value::Vector{Float64})
    if has_core_point(normalization)
        check_reverse_polar_dimension(normalization.core_point_x, length(x_value), "core_point_x")
        check_reverse_polar_dimension(normalization.core_point_t, length(t_value), "core_point_t")
        return x_value .- normalization.core_point_x, t_value .- normalization.core_point_t
    end

    if has_core_direction(normalization)
        is_default_vertical_direction(normalization) &&
            return zeros(length(x_value)), ones(length(t_value))

        check_reverse_polar_dimension(normalization.core_diretion_x, length(x_value), "core_diretion_x")
        check_reverse_polar_dimension(normalization.core_diretion_t, length(t_value), "core_diretion_t")
        return normalization.core_diretion_x, normalization.core_diretion_t
    end

    throw(ArgumentError("ReversePolarNormalization requires a core point, a core direction, or the default vertical direction."))
end

function dcglp_reverse_polar_direction(normalization::ReversePolarNormalization, x_value::Vector{Float64}, t_value::Vector{Float64})
    direction_x, direction_t = reverse_polar_direction(normalization, x_value, t_value)
    has_core_point(normalization) && return direction_x, direction_t
    return .-direction_x, .-direction_t
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
end

function should_fallback_typical(normalization::ReversePolarNormalization, oracle::SplitOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    has_core_point(normalization) || return false

    direction_x, direction_t = reverse_polar_direction(normalization, x_value, t_value)

    if LinearAlgebra.norm(vcat(direction_x, direction_t), Inf) <= oracle.param.zero_tol
        return true
    end

    return false
end

function update_dcglp_for_candidate!(normalization::ReversePolarNormalization, oracle::SplitOracle, x_value::Vector{Float64}, t_value::Vector{Float64})
    set_normalized_rhs.(oracle.dcglp[:conx], x_value)
    set_normalized_rhs.(oracle.dcglp[:cont], t_value)

    direction_x, direction_t = dcglp_reverse_polar_direction(normalization, x_value, t_value)
    update_reverse_polar_constraints!(oracle, direction_x, direction_t)
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
    direction_x, direction_t = dcglp_reverse_polar_direction(normalization, x_value, reference_t)
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
end

function disjunctive_cut_normalization_value(
    normalization::ReversePolarNormalization,
    gamma_x::Vector{Float64},
    gamma_t::Vector{Float64},
    common::SplitOracleParam,
    current_lb::Float64,
    x_value::Vector{Float64},
    t_value::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int},
)
    direction_x, direction_t = dcglp_reverse_polar_direction(normalization, x_value, t_value)
    direction_value = dot(gamma_x, direction_x) + dot(gamma_t, direction_t)
    abs(direction_value) > common.zero_tol ||
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
    has_core_point(normalization) ||
        throw(ArgumentError("`set_core_point!` requires a directional ReversePolarNormalization."))
    length(core_point_x) == length(normalization.core_point_x) ||
        throw(DimensionMismatch("`core_point_x` has length $(length(core_point_x)) but expected $(length(normalization.core_point_x))."))
    length(core_point_t) == length(normalization.core_point_t) ||
        throw(DimensionMismatch("`core_point_t` has length $(length(core_point_t)) but expected $(length(normalization.core_point_t))."))

    normalization.core_point_x .= core_point_x
    normalization.core_point_t .= core_point_t
end

function update_reverse_polar_constraints!(
    oracle::SplitOracle,
    direction_x::Vector{Float64},
    direction_t::Vector{Float64},
)
    dcglp = oracle.dcglp

    if haskey(dcglp, :con_reverse_polar_x) && haskey(dcglp, :con_reverse_polar_t) &&
       length(dcglp[:con_reverse_polar_x]) == length(direction_x) && length(dcglp[:con_reverse_polar_t]) == length(direction_t)
        tau = dcglp[:tau]
        for j in eachindex(direction_x)
            set_normalized_coefficient(dcglp[:con_reverse_polar_x][j], tau, direction_x[j])
        end
        for j in eachindex(direction_t)
            set_normalized_coefficient(dcglp[:con_reverse_polar_t][j], tau, direction_t[j])
        end
        return nothing
    end

    delete_registered_constraints!(dcglp, :con_reverse_polar_x)
    delete_registered_constraints!(dcglp, :con_reverse_polar_t)

    @constraint(
        dcglp,
        con_reverse_polar_x[j in eachindex(direction_x)],
        dcglp[:sx][j] + direction_x[j] * dcglp[:tau] == 0.0,
    )
    @constraint(
        dcglp,
        con_reverse_polar_t[j in eachindex(direction_t)],
        dcglp[:st][j] + direction_t[j] * dcglp[:tau] == 0.0,
    )
end
