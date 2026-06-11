"""
    apply_lift_or_strengthen(dcglp, gamma_x, zero_indices, one_indices; lift, strengthen, zero_tol, gamma_0=0.0)

Apply lifting (preferred if `lift` and any zero/one index exists) or pure
strengthening to the cut coefficients. Returns the (possibly modified)
`(gamma_x, gamma_0)` pair.
"""
function apply_lift_or_strengthen(
    dcglp::Model,
    gamma_x::Vector{Float64},
    zero_indices::Vector{Int},
    one_indices::Vector{Int};
    lift::Bool,
    strengthen::Bool,
    zero_tol::Float64,
    gamma_0::Float64 = 0.0,
)
    if lift && (!isempty(zero_indices) || !isempty(one_indices))
        return lift_x_coefficients(
            dcglp, gamma_x, gamma_0, zero_indices, one_indices;
            strengthen = strengthen, zero_tol = zero_tol,
        )
    elseif strengthen
        sigma, delta = read_strengthening_duals(dcglp)
        return -strengthen_coefficients(-gamma_x, sigma, delta; zero_tol = zero_tol), gamma_0
    end
    return gamma_x, gamma_0
end

read_strengthening_duals(dcglp::Model) = (
    Dict(1 => dual(dcglp[:con_split_kappa]), 2 => dual(dcglp[:con_split_nu])),
    Dict(1 => dual.(dcglp[:condelta][1, :]), 2 => dual.(dcglp[:condelta][2, :])),
)

"""
    lift_x_coefficients(dcglp, gamma_x, gamma_0, zero_indices, one_indices; strengthen, zero_tol)

Apply lifting (and optional strengthening) to `gamma_x`. Returns the lifted
`(gamma_x, gamma_0)` pair in the common cut orientation.
"""
function lift_x_coefficients(
    dcglp::Model,
    gamma_x::Vector{Float64},
    gamma_0::Float64,
    zero_indices::Vector{Int},
    one_indices::Vector{Int};
    strengthen::Bool,
    zero_tol::Float64,
)
    zeta_k = !isempty(zero_indices) ? dual.(dcglp[:con_zeta][1, :]) : Float64[]
    zeta_v = !isempty(zero_indices) ? dual.(dcglp[:con_zeta][2, :]) : Float64[]
    xi_k = !isempty(one_indices) ? dual.(dcglp[:con_xi][1, :]) : Float64[]
    xi_v = !isempty(one_indices) ? dual.(dcglp[:con_xi][2, :]) : Float64[]

    lifted_gamma_0 = gamma_0 - sum(max.(xi_k, xi_v))
    lifted_gamma_x = -gamma_x
    lifted_gamma_x[zero_indices] .= -gamma_x[zero_indices] .+ max.(zeta_k, zeta_v)
    lifted_gamma_x[one_indices] .= -gamma_x[one_indices] .- max.(xi_k, xi_v)

    if strengthen
        sigma, base_delta = read_strengthening_duals(dcglp)
        delta_1 = copy(base_delta[1])
        delta_2 = copy(base_delta[2])
        delta_1[zero_indices] .+= -zeta_k .+ max.(zeta_k, zeta_v)
        delta_2[zero_indices] .+= -zeta_v .+ max.(zeta_k, zeta_v)
        lifted_gamma_x = strengthen_coefficients(lifted_gamma_x, sigma, Dict(1 => delta_1, 2 => delta_2); zero_tol = zero_tol)
    end

    return -lifted_gamma_x, lifted_gamma_0
end

function normalize_directional_duals!(
    gamma_x::AbstractVector{Float64},
    gamma_t::AbstractVector{Float64},
    direction_x::Vector{Float64},
    direction_t::Vector{Float64};
    zero_tol::Float64,
)
    direction_value = dot(gamma_x, direction_x) + dot(gamma_t, direction_t)
    abs(direction_value) > zero_tol ||
        throw(AlgorithmException("DirectionalReversePolarNormalization cut normalization failed because the directional support is numerically zero."))
    gamma_x ./= direction_value
    gamma_t ./= direction_value
    return direction_value
end

function strengthen_coefficients(gamma_x, sigma, delta; zero_tol = 1.0e-9)
    a_1 = gamma_x .- delta[1]
    a_2 = gamma_x .- delta[2]
    sigma_sum = sigma[1] + sigma[2]
    sigma_sum < zero_tol && return copy(gamma_x)
    m = (a_1 .- a_2) ./ sigma_sum
    return min.(a_1 .- sigma[1] .* floor.(m), a_2 .+ sigma[2] .* ceil.(m))
end

function compute_norm_value(gamma_x, gamma_t, norm::AbstractNorm)
    if norm.p == 1.0
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), Inf)
    elseif norm.p == 2.0
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), 2.0)
    elseif norm.p == Inf
        norm_value = LinearAlgebra.norm(vcat(gamma_x, gamma_t), 1.0)
    else
        throw(UndefError("Unsupported norm type: $(typeof(norm))"))
    end
    return max(1.0, norm_value)
end
