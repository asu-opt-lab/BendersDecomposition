# Single iteration information
"""
    DCGLPIterationInfo

Stores information about a single iteration of the DCGLP (Disjunctive Cut Generating Linear Program) algorithm.

# Fields
- `iter::Int`: Current iteration number
- `LB::Float64`: Lower bound
- `UB::Float64`: Upper bound
- `UB_k::Float64`: Upper bound for k-subproblem
- `UB_v::Float64`: Upper bound for v-subproblem
- `gap::Float64`: Optimality gap
- `master_time::Float64`: Time spent in master problem
- `sub_k_time::Float64`: Time spent in k-subproblem
- `sub_v_time::Float64`: Time spent in v-subproblem
- `total_time::Float64`: Total elapsed time
"""
struct DCGLPIterationInfo
    iter::Int
    LB::Float64
    UB::Float64
    UB_k::Float64
    UB_v::Float64
    gap::Float64
    master_time::Float64
    sub_k_time::Float64
    sub_v_time::Float64
    total_time::Float64
end

# Stores the current state of the Benders algorithm
"""
    DCGLPState

Maintains the current state of the Benders algorithm for DCGLP.

# Fields
- `iteration::Int`: Current iteration count
- `LB::Float64`: Current lower bound
- `UB::Float64`: Current upper bound
- `UB_k::Float64`: Current upper bound for k-subproblem
- `UB_v::Float64`: Current upper bound for v-subproblem
- `gap::Float64`: Current optimality gap
"""
mutable struct DCGLPState
    iteration::Int
    LB::Float64
    UB::Float64
    UB_k::Float64
    UB_v::Float64
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

"""
    DCGLPIterationLog

Tracks the history and timing information for the DCGLP algorithm.

# Fields
- `iterations::Vector{DCGLPIterationInfo}`: History of all iterations
- `start_time::Float64`: Algorithm start time
- `master_time::Float64`: Cumulative time spent in master problem
- `sub_k_time::Float64`: Cumulative time spent in k-subproblem
- `sub_v_time::Float64`: Cumulative time spent in v-subproblem
"""
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



function generate_strengthened_cuts(dcglp::DCGLP, ::PureDisjunctiveCut)
    optimize!(dcglp.model)
    γₜ = dual(dcglp.model[:cont])
    γ₀ = dual(dcglp.model[:con0])
    γₓ = dual.(dcglp.model[:conx])
    return γ₀, γₓ, γₜ
end

"""
    generate_strengthened_cuts(dcglp::DCGLP, ::StrengthenedDisjunctiveCut)

Generate strengthened disjunctive cuts.

Returns a tuple (γ₀, γₓ, γₜ) representing coefficients for the constant term,
decision variables, and time variable respectively. This method computes stronger
valid inequalities using dual information.
"""
function generate_strengthened_cuts(dcglp::DCGLP, ::StrengthenedDisjunctiveCut)
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





