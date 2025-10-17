export DualDecompositionParam, DualDecompositionLog, DualDecomposition

include("$(dirname(dirname(@__DIR__)))/example/cflp/oracle.jl")

mutable struct DualDecompositionParam <: AbstractOracleParam
    # Uniform distribution part
    max_iter::Float64
    obj_limit::Float64
    halving_value::Float64
    step_size::Float64
    γ_max::Float64

    # Diminihisng part
    stepsize_bound::Float64
    stepsize_constant::Float64

    # for uniform distribution
    function DualDecompositionParam(; max_iter = 100, obj_limit = +Inf, stepsize = 0.1, halving_value = 1.0, γ_max = max_iter*0.2, stepsize_bound = 1e-3, stepsize_constant = 0.1)
        new(max_iter, obj_limit, stepsize, halving_value, γ_max, stepsize_bound, stepsize_constant)
    end
    # for diminishing
    # function DualDecompositionParam(; max_iter = 100, obj_limit = +Inf, stepsize = 0.1, halving_value = 1.0, γ_max = max_iter*0.2, stepsize_bound = 1e-2, stepsize_constant = 0.04)
    #     # for α = C/√k, use stepsize constant as 0.1
    #     new(max_iter, obj_limit, stepsize, halving_value, γ_max, stepsize_bound, stepsize_constant)
    # end
end

mutable struct DualDecompositionLog <: AbstractOracleLog
    pri_var::Matrix{Float64}
    dual_var::Dict{Symbol, Array{Float64}}
    α_set::Vector{Any}
    λ_k_diff_set::Vector{Any}
    prev_is_in_L::Bool
    dual_obj::Float64

    # Uniform distribution part
    vogel_time::Float64

    function DualDecompositionLog(data::Data; α_set = [], λ_k_diff_set = [], prev_is_in_L = false, dual_obj = -Inf)
        m = data.dim_x; n = data.problem.n_customers
        pri_var =  zeros(Float64, m, n); dual_var = Dict(:λ => zeros(m), :σ => zeros(n), :δ => zeros(Float64, m, n))
        new(pri_var, dual_var, α_set, λ_k_diff_set, prev_is_in_L, dual_obj)
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

# update order: y → λ (uniform distribution)
function generate_cuts(oracle::DualDecomposition, x::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    println("DualDecomposition")
    data = oracle.data; param = oracle.oracle_param; log = oracle.oracle_log; typical_oracle = oracle.typical_oracle

    log.vogel_time = @elapsed param.obj_limit = UB_approximation(data, x)

    # Just for comparison with optimal info
    _, _, f_x = generate_cuts(typical_oracle, x, t_value)
    opt_y = value.(typical_oracle.model[:y]); opt_λ = dual.(typical_oracle.model[:capacity]); @debug opt_λ, f_x

    # call necessary parameters
    y_k =log.pri_var; λ_k = log.dual_var[:λ]; halving_value = copy(param.halving_value); α = copy(param.step_size)
    d = data.problem.demands; u= data.problem.capacities; c = data.problem.costs .* d'

    iter, γ = 1, 0; y_time =0; res_time = 0; λ_time = 0; y_dual_time = 0; residual = Vector{Float64}(undef, length(u))
    while α >= param.stepsize_bound
        # Update y, compute residual and LB
        s1 = time()
        y_k, _, _, _, _ = update_y(y_k, λ_k, c, d, x)
        y_time = (time() -s1)

        s2 = time()
        @threads for i in eachindex(u)
            residual[i] = dot(d, @view y_k[i, :]) - u[i] * x[i]
        end
        z_lb = sum(c .* y_k) + λ_k' * residual
        res_time = (time() - s2)

        # Update step size
        α = halving_value * abs(param.obj_limit - z_lb)/norm(residual)^2 # Polyak's stepsize

        # Update dual value
        s3 = time()
        λ_k = update_λ(λ_k, y_k, u, d, x, α)
        λ_time = (time() - s3)

        # Update halving parameter
        γ += 1
        γ > param.γ_max && (halving_value = halving_value * 0.5; γ = 0)

        iter += 1
    end
    @info iter
    s4 = time()
    y_k, sorted_indices, c_sorted, _, critical_facility_indices = update_y(y_k, λ_k, c, d, x)
    retrieve_dual_values(log, d, sorted_indices, c_sorted, critical_facility_indices)
    y_dual_time = (time() - s4)
    println("y_time: $y_time, res_time: $res_time, λ_time: $λ_time, y_dual_time: $y_dual_time")
    
    feasibility = ((log.dual_var[:σ]') .- log.dual_var[:δ] .- (log.dual_var[:λ] * d')) .<= c .+ 1e-9
    @assert all(feasibility) "Dual feasibility violated"

    dual_obj = sum(log.dual_var[:σ]) - dot(x, vec(sum(log.dual_var[:δ],dims = 2))) - dot(log.dual_var[:λ], u.*x)

    @assert dual_obj <= f_x[1] + 1e-3 "dual($dual_obj) > f^*($(f_x[1]))"
    log.dual_obj = dual_obj
    if dual_obj >= t_value[1]* (1 + 1e-9) + 1e-6/tol_normalize
        a_x = -log.dual_var[:λ].*u .- [sum(log.dual_var[:δ][i,:]) for i in eachindex(u)]
        a_t = [-1.0]
        a_0 = sum(log.dual_var[:σ]) 
        # return false, [Hyperplane(a_x, a_t, a_0)], [NaN]
        return false, [Hyperplane(a_x, a_t, a_0)], f_x
    else
        typical_oracle.oracle_param.dual_value[:λ] = -λ_k; typical_oracle.oracle_param.dual_value[:δ] = log.dual_var[:δ]; typical_oracle.oracle_param.dual_value[:σ] = log.dual_var[:σ]
        is_in_L, hyperplanes, f_x = generate_cuts(typical_oracle, x, t_value; time_limit = time_limit)
        println("DD_dual_obj ≤ t")
        return is_in_L, hyperplanes, f_x
    end
end

# update order: y → λ (diminishing)
# function generate_cuts(oracle::DualDecomposition, x::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
#     data = oracle.data; param = oracle.oracle_param; log = oracle.oracle_log; typical_oracle = oracle.typical_oracle

#     # call necessary parameters
#     y_k =log.pri_var; λ_k = log.dual_var[:λ]
#     d = data.problem.demands; u= data.problem.capacities; c = data.problem.costs .* d'

#     iter = 1; y_time =0; λ_time = 0; y_dual_time = 0
#     while true
#         # Update y, compute residual and LB
#         s1 = time()
#         y_k, _, _, _, _ = update_y(y_k, λ_k, c, d, x)
#         y_time += (time() - s1)

#         # Update step size
#         # α = param.stepsize_constant/sqrt(iter)
#         α = param.stepsize_constant/iter^0.25

#         # Update dual value
#         s2 = time()
#         λ_k = update_λ(λ_k, y_k, u, d, x, α)
#         λ_time += (time() - s2)

#         # Store necessary results
#         # push!(log.α_set, α); # push!(log.λ_k_diff_set, norm(opt_λ .+ λ_k))

#         # Termination
#         α < param.stepsize_bound && break
#         iter +=1
#     end
#     @info iter
#     s3 = time()
#     y_k, sorted_indices, c_sorted, _, critical_facility_indices = update_y(y_k, λ_k, c, d, x)
#     retrieve_dual_values(log, d, sorted_indices, c_sorted, critical_facility_indices)
#     y_dual_time += (time() - s3)
    
#     feasibility = ((log.dual_var[:σ]') .- log.dual_var[:δ] .- (log.dual_var[:λ] * d')) .<= c .+ 1e-9
#     @assert all(feasibility) "Dual feasibility violated"
    
#     dual_obj = sum(log.dual_var[:σ]) - dot(x, vec(sum(log.dual_var[:δ],dims = 2))) - dot(log.dual_var[:λ], u.*x)

#     println("y_time: $y_time, λ_time: $λ_time, y_dual_time: $y_dual_time")
#     log.dual_obj = dual_obj
#     if dual_obj >= t_value[1]* (1 + 1e-6) + 1e-6/tol_normalize
#         a_x = -log.dual_var[:λ].*u .- [sum(log.dual_var[:δ][i,:]) for i in eachindex(u)]
#         a_t = [-1.0]
#         a_0 = sum(log.dual_var[:σ]) 
#         return false, [Hyperplane(a_x, a_t, a_0)], [NaN]
#     else
#         is_in_L, hyperplanes, f_x = generate_cuts(typical_oracle, x, t_value; time_limit = time_limit)
#         println("DD_dual_obj ≤ t")
#         return is_in_L, hyperplanes, f_x
#     end
# end