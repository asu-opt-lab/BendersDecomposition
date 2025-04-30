
"""
prototype for cutting-plane algorithms for DCGLP
"""
function solve_dcglp!(oracle::AbstractDisjunctiveOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; rtol = 1e-8, atol = 1e-9, time_limit = time_limit, throw_typical_cuts_for_errors = true, include_disjuctive_cuts_to_hyperplanes = true)
    throw(UndefError("Update `solve_dcglp!` for $(typeof(oracle))"))
end

function generate_disjunctive_cut(dcglp::Model; strengthen = false, zero_tol = 1e-7)
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    gamma_0 = dual(dcglp[:con0])

    if strengthen 
        sigma = Dict(1 => dual(dcglp[:con_split_kappa]), 2 => dual(dcglp[:con_split_nu]))
        delta = Dict(1 => dual.(dcglp[:condelta][1,:]), 2 => dual.(dcglp[:condelta][2,:]))
        gamma_x = -strengthening!(-gamma_x, sigma, delta; zero_tol = zero_tol)
    end
    
    return gamma_x, gamma_t, gamma_0
end

function generate_lifted_disjunctive_cut(dcglp::Model, norm::LpNorm, zero_indices::Vector{Int64}, one_indices::Vector{Int64}; strengthen = false, zero_tol = 1e-7)
    gamma_x = dual.(dcglp[:conx])
    gamma_t = dual.(dcglp[:cont])
    gamma_0 = dual(dcglp[:con0])

    zeta_k = !isempty(zero_indices) ? dual.(dcglp[:con_zeta][1,:]) : Float64[]
    zeta_v = !isempty(zero_indices) ? dual.(dcglp[:con_zeta][2,:]) : Float64[] 
    xi_k = !isempty(one_indices) ? dual.(dcglp[:con_xi][1,:]) : Float64[] 
    xi_v = !isempty(one_indices) ? dual.(dcglp[:con_xi][2,:]) : Float64[] 

    # coefficients for lifted cut
    lifted_gamma_0 = gamma_0 - sum(max.(xi_k, xi_v))
    lifted_gamma_x = zeros(Float64, length(gamma_x))
    lifted_gamma_x .= -gamma_x

    lifted_gamma_x[zero_indices] = -gamma_x[zero_indices] .+ max.(zeta_k, zeta_v)
    lifted_gamma_x[one_indices] = -gamma_x[one_indices] .- max.(xi_k, xi_v)

    if strengthen
        sigma = Dict(1 => dual(dcglp[:con_split_kappa]), 2 => dual(dcglp[:con_split_nu]))
        delta_1 = dual.(dcglp[:condelta][1,:])
        delta_2 = dual.(dcglp[:condelta][2,:])
        delta_1[zero_indices] += (-zeta_k + max.(zeta_k, zeta_v))
        delta_2[zero_indices] += (-zeta_v + max.(zeta_k, zeta_v)) 
        delta = Dict(1 => delta_1, 2 => delta_2)
        lifted_gamma_x = strengthening!(lifted_gamma_x, sigma, delta; zero_tol = zero_tol)
    end

    # compute normalization value
    norm_value = compute_norm_value(lifted_gamma_x, gamma_t, norm)

    return (-lifted_gamma_x, gamma_t, lifted_gamma_0) ./ norm_value
end

function strengthening!(gamma_x, sigma, delta; zero_tol = 1e-7)
    @debug "dcglp strengthening - sigma values: [σ₁: $(sigma[1]), σ₂: $(sigma[2])]"
    @debug "dcglp strengthening - delta values: [δ₁: $(delta[1]), δ₂: $(delta[2])]"
    
    a₁ = gamma_x .- delta[1]
    a₂ = gamma_x .- delta[2]
    sigma_sum = sigma[1] + sigma[2]
    if sigma_sum >= zero_tol
        m = (a₁ .- a₂) / sigma_sum
        m_lb = floor.(m)
        m_ub = ceil.(m)
        gamma_x = min.(a₁-sigma[1]*m_lb, a₂+sigma[2]*m_ub)
    end
    return gamma_x
end


include("DcglpGammaNorm.jl")
include("DcglpStdNorm.jl")
