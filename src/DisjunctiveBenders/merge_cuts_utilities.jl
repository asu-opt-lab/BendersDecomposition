function generate_disjunctive_cuts(dcglp::DCGLP, ::PureDisjunctiveCut)

    γₜ = -dual.(dcglp.model[:cont])
    γ₀ = dual(dcglp.model[:con0])
    γₓ = -dual.(dcglp.model[:conx])
    return γ₀, γₓ, γₜ
end


function generate_disjunctive_cuts(dcglp::DCGLP, ::StrengthenedDisjunctiveCut)
    
    σ₁::Float64 = dual(dcglp.disjunctive_inequalities_constraints[1])
    σ₂::Float64 = dual(dcglp.disjunctive_inequalities_constraints[2])
    γₜ::Float64 = dual(dcglp.γ_constraints[:γₜ])
    γ₀::Float64 = dual(dcglp.γ_constraints[:γ₀])
    γₓ::Vector{Float64} = -dual.(dcglp.γ_constraints[:γₓ])

    if abs(σ₁) <= 1e-6 && abs(σ₂) <= 1e-6
        return γ₀, -γₓ, γₜ
    end

    γ₁ = γₓ .- dual.(dcglp.model[:conv1])
    γ₂ = γₓ .- dual.(dcglp.model[:conv2])
    σ_sum = σ₂ + σ₁
    if !iszero(σ_sum)
        m = (γ₁ .- γ₂) / σ_sum
        m_lb = floor.(m)
        m_ub = ceil.(m)
        γₓ = min.(γ₁-σ₁*m_lb, γ₂+σ₂*m_ub)
    end
    
    return γ₀, -γₓ, γₜ
end