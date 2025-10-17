export DualDecompositionParam, DualDecompositionLog, DualDecomposition

include("$(dirname(dirname(@__DIR__)))/example/cflp/oracle.jl")

mutable struct DualDecompositionParam <: AbstractOracleParam
    stale_lim::Int64
    stepsize_bound::Float64
    stepsize_constant::Float64

    function DualDecompositionParam(; stale_lim = 3, stepsize_bound = 1e-2, stepsize_constant = 0.1)
        new(stale_lim, stepsize_bound, stepsize_constant)
    end
end

mutable struct DualDecompositionLog <: AbstractOracleLog
    pri_var::Matrix{Float64}
    dual_var::Dict{Symbol, Array{Float64}}
    α_set::Vector{Any}
    λ_k_diff_set::Vector{Any}
    prev_is_in_L::Bool
    fin_LB_set::Vector{Any}

    function DualDecompositionLog(data::Data; α_set = [], λ_k_diff_set = [], prev_is_in_L = false, fin_LB_set = [])
        m = data.dim_x; n = data.problem.n_customers
        pri_var =  zeros(Float64, m, n); dual_var = Dict(:λ => zeros(m), :σ => zeros(n), :δ => zeros(Float64, m, n))
        new(pri_var, dual_var, α_set, λ_k_diff_set, prev_is_in_L, fin_LB_set)
    end
end

mutable struct DualDecomposition <: AbstractTypicalOracle
    data::Data
    typical_oracle::AbstractOracle
    oracle_param::DualDecompositionParam
    oracle_log::DualDecompositionLog
    flag_bnb::Bool

    function DualDecomposition(data::Data; typical_oracle::AbstractOracle = ClassicalOracle(data),
        oracle_param::DualDecompositionParam = DualDecompositionParam(), oracle_log::DualDecompositionLog = DualDecompositionLog(data), flag_bnb = false)
        @debug "Building Dual Decomposition Oracle"
        new(data, typical_oracle, oracle_param, oracle_log, flag_bnb)
    end

    DualDecomposition() = new()
end

# # update order: y → λ
function generate_cuts(oracle::DualDecomposition, x::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    data = oracle.data; param = oracle.oracle_param; log = oracle.oracle_log; typical_oracle = oracle.typical_oracle

    if log.prev_is_in_L || (length(log.fin_LB_set) >= param.stale_lim + 1 && all(isapprox.(log.fin_LB_set[end], log.fin_LB_set[end-param.stale_lim:end-1])))
        println("Exact solver")
        log.prev_is_in_L = false
        is_in_L, hyperplanes, f_x = generate_cuts(typical_oracle, x, t_value; time_limit = time_limit)
        return is_in_L, hyperplanes, f_x
    else
        println("Dual Decomposition")

        # call necessary parameters
        y_k =log.pri_var; λ_k = log.dual_var[:λ]
        d = data.problem.demands; u= data.problem.capacities; c = data.problem.costs .* d'

        iter = 1
        while true
            # Update y, compute residual and LB
            s1 = time()
            y_k, _, _, _, _ = update_y(y_k, λ_k, c, d, x)
            y_time += (time() - s1)

            # Update step size
            α = param.stepsize_constant/sqrt(iter)

            # Update dual value
            s2 = time()
            λ_k = update_λ(λ_k, y_k, u, d, x, α)
            λ_time += (time() - s2)

            # Store necessary results
            # push!(log.α_set, α); # push!(log.λ_k_diff_set, norm(opt_λ .+ λ_k))

            # Termination
            α < param.stepsize_bound && break
            iter +=1
        end
        s3 = time()
        y_k, sorted_indices, c_sorted, _, critical_facility_indices = update_y(y_k, λ_k, c, d, x)
        retrieve_dual_values(log, d, sorted_indices, c_sorted, critical_facility_indices)
        y_dual_time += (time() - s3)
        
        feasibility = ((log.dual_var[:σ]') .- log.dual_var[:δ] .- (log.dual_var[:λ] * d')) .<= c .+ 1e-9
        @assert all(feasibility) "Dual feasibility violated"
        
        dual_obj = sum(log.dual_var[:σ]) - dot(x, vec(sum(log.dual_var[:δ],dims = 2))) - dot(log.dual_var[:λ], u.*x)
        !oracle.flag_bnb && push!(log.fin_LB_set, dual_obj)

        println("y_time: $y_time, λ_time: $λ_time, y_dual_time: $y_dual_time")
        log.prev_is_in_L = false
        if dual_obj >= t_value[1]* (1 + 1e-6) + 1e-6/tol_normalize
            a_x = -log.dual_var[:λ].*u .- [sum(log.dual_var[:δ][i,:]) for i in eachindex(u)]
            a_t = [-1.0]
            a_0 = sum(log.dual_var[:σ]) 
            return false, [Hyperplane(a_x, a_t, a_0)], [NaN]
        else
            if !oracle.flag_bnb
                log.prev_is_in_L = true
                return false, [Hyperplane(zeros(length(x)), [0.0], 0.0)], [NaN]
            else
                is_in_L, hyperplanes, f_x = generate_cuts(typical_oracle, x, t_value; time_limit = time_limit)
                println("lambda_DD:$(log.dual_var[:λ])")
                return is_in_L, hyperplanes, f_x
            end
        end
    end
end