
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
    LB_set::Vector{Float64}

    # Constructor with default values
    function DCGLPState()
        new(0, -Inf, Inf, Inf, Inf, Inf, [])
    end

    # Constructor with specified values
    function DCGLPState(LB::Float64, UB::Float64, UB_k::Float64, UB_v::Float64, gap::Float64, iteration::Int = 0)
        new(iteration, UB, LB, UB_k, UB_v, gap, [])
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








