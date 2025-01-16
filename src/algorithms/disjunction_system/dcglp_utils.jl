
struct DCGLPIterationInfo
    iter::Int
    LB::Float64
    UB::Float64
    UB_k::Union{Vector{Float64}, Float64}
    UB_v::Union{Vector{Float64}, Float64}
    gap::Float64
    master_time::Float64
    sub_k_time::Float64
    sub_v_time::Float64
    total_time::Float64
end


mutable struct DCGLPState
    iteration::Int
    LB::Float64
    UB::Float64
    UB_k::Union{Vector{Float64}, Float64}
    UB_v::Union{Vector{Float64}, Float64}
    gap::Float64

    # Constructor with default values
    function DCGLPState()
        new(0, -Inf, Inf, Inf, Inf, Inf)
    end

    # Constructor with specified values
    function DCGLPState(LB::Float64, UB::Float64, UB_k::Float64, UB_v::Float64, gap::Float64, iteration::Int = 0)
        new(iteration, UB, LB, UB_k, UB_v, gap)
    end
end


mutable struct DCGLPIterationLog
    iterations::Vector{DCGLPIterationInfo}
    start_time::Float64
    master_time::Float64
    sub_k_time::Float64
    sub_v_time::Float64

    function DCGLPIterationLog()
        new(DCGLPIterationInfo[], time(), 0.0, 0.0, 0.0)
    end
end


function get_total_time(log::DCGLPIterationLog)
    return time() - log.start_time
end


function record_iteration!(log::DCGLPIterationLog, state::DCGLPState)
    push!(log.iterations, DCGLPIterationInfo(
        state.iteration,
        state.LB,
        state.UB,
        state.UB_k,
        state.UB_v,
        state.gap,
        log.master_time,
        log.sub_k_time,
        log.sub_v_time,
        get_total_time(log)
    ))
end




# ============================================================================
# Generate strengthened cuts for the DCGLP
# ============================================================================



function generate_disjunctive_cuts(dcglp::DCGLP, ::PureDisjunctiveCut)
    optimize!(dcglp.model)
    γₜ = dual.(dcglp.model[:cont])
    γ₀ = dual(dcglp.model[:con0])
    γₓ = dual.(dcglp.model[:conx])
    return γ₀, γₓ, γₜ
end


function generate_disjunctive_cuts(dcglp::DCGLP, ::StrengthenedDisjunctiveCut)
    optimize!(dcglp.model)
    
    σ₁::Float64 = dual(dcglp.disjunctive_inequalities_constraints[1])
    σ₂::Float64 = dual(dcglp.disjunctive_inequalities_constraints[2])
    γₜ::Float64 = dual(dcglp.γ_constraints[:γₜ])
    γ₀::Float64 = dual(dcglp.γ_constraints[:γ₀])
    γₓ::Vector{Float64} = dual.(dcglp.γ_constraints[:γₓ])
    if iszero(σ₁) && iszero(σ₂)
        return γ₀, γₓ, γₜ
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
    
    return γ₀, γₓ, γₜ
end





