export DualDecompositionParam, DualDecompositionLog, DualDecomposition

mutable struct DualDecompositionParam <: AbstractOracleParam
    LB_stag_consecutive_iter::Float64
    max_iter::Float64
    obj_limit::Float64
    step_size::Float64
    halving_value::Float64
    γ_max::Float64
    stepsize_bound::Float64

    function DualDecompositionParam(; LB_stag_consecutive_iter = 3, max_iter = 1000, obj_limit = +Inf, step_size = 0.1, halving_value = 2.0, γ_max = max_iter*0.1, stepsize_bound = 1e-4)
        new(LB_stag_consecutive_iter, max_iter, obj_limit, step_size, halving_value, γ_max, stepsize_bound)
    end
end

mutable struct DualDecompositionLog <: AbstractOracleLog
    pri_var::Matrix{Float64}
    slack::Vector{Float64}
    dual_var::Dict{Symbol, Array{Float64}}
    best_obj::Float64
    obj_set::Vector{Any}
    α_set::Vector{Any}
    λ_k_norm_set::Vector{Any}
    λ_k_diff_set::Vector{Any}
    residual_norm_set::Vector{Any}
    approx_obj_set::Vector{Any}
    flag_exceed::Bool

    function DualDecompositionLog(data::Data; best_obj = -Inf, obj_set = [], α_set = [], λ_k_norm_set = [], λ_k_diff_set = [], residual_norm_set = [], approx_obj_set = [], flag_exceed = false)
        pri_var =  zeros(Float64, data.dim_x, data.problem.n_customers); slack = zeros(data.dim_x)
        dual_var = Dict(:λ => zeros(data.dim_x), :σ => zeros(data.problem.n_customers), :δ => zeros(Float64, data.dim_x, data.problem.n_customers))
        new(pri_var, slack, dual_var, best_obj, obj_set, α_set, λ_k_norm_set, λ_k_diff_set, residual_norm_set, approx_obj_set, flag_exceed)
    end
end

mutable struct DualDecomposition <: AbstractPrimalDualOracle
    
    oracle_param::DualDecompositionParam
    oracle_log::DualDecompositionLog

    function DualDecomposition(data::Data; oracle_param::DualDecompositionParam = DualDecompositionParam(), oracle_log::DualDecompositionLog = DualDecompositionLog(data))
        @debug "Building Dual Decomposition Oracle"
        new(oracle_param, oracle_log)
    end

    DualDecomposition() = new()
end

# update y_k first, then update lambda
# function generate_cuts(data::Data, oracle::DualDecomposition, x::Vector{Float64}, f_x, opt_λ, global_iter)
#     param = oracle.oracle_param; log = oracle.oracle_log
#     y_k =log.pri_var; λ_k = log.dual_var[:λ]; best_obj = copy(log.best_obj); halving_value = copy(param.halving_value); α = copy(param.step_size)
#     d, u = data.problem.demands, data.problem.capacities; c = data.problem.costs .* d'

#     iter, γ = 1, 0
#     # while α >= param.stepsize_bound
#     while α >= 1e-2
#         # Update y, compute residual and LB
#         y_k, _, _, _, _ = update_y(y_k, λ_k, c, d, x)
#         residual = [d'*y_k[i,:] - u[i]*x[i] for i in eachindex(u)]
#         z_lb = sum(c .* y_k) + λ_k' * residual

#         # Update step size
#         # α = halving_value * abs(param.obj_limit - z_lb)/norm(residual)^2 # Polyak's stepsize
#         α = 0.1/sqrt(iter) # diminishing

#         # Update dual value
#         λ_k = update_λ(λ_k, y_k, u, d, x, α)

#         # Keep the best solution and obj
#         if best_obj < z_lb 
#             best_obj = z_lb
#         else
#             γ += 1
#             if γ > param.γ_max
#                 halving_value = halving_value * 0.5
#                 γ = 0
#             end
#         end
#         push!(log.α_set, α); push!(log.obj_set, norm(f_x - z_lb)); push!(log.λ_k_diff_set, norm(opt_λ .+ λ_k))
#         # push!(log.λ_k_norm_set, norm(λ_k)); push!(log.residual_norm_set, norm(residual)); push!(log.approx_obj_set, sum(c .* y_k))

#         # Handle divergence of Polyak caused by not diminishing stepsize (do not immediately terminate)
#         (best_obj > param.obj_limit && iter > param.max_iter) && (@info "best f: $(best_obj) > overestimated f: $(param.obj_limit)"; log.flag_exceed = true; break)
#         iter += 1
#     end
#     _, sorted_indices, c_sorted, _, critical_facility_indices = update_y(y_k, λ_k, c, d, x)
#     retrieve_dual_values(log, d, sorted_indices, c_sorted, critical_facility_indices)

#     a_x = -log.dual_var[:λ].*u .- [sum(log.dual_var[:δ][i,:]) for i in eachindex(u)]
#     a_t = [-1.0]
#     a_0 = sum(log.dual_var[:σ]) 

#     # gen_figure(log, global_iter)
#     # empty!(log.α_set); empty!(log.obj_set); empty!(log.λ_k_diff_set)

#     @assert (sum(log.dual_var[:σ]) - dot(x, vec(sum(log.dual_var[:δ],dims = 2))) - dot(log.dual_var[:λ], u.*x) <= f_x +1e-3) "Dual obj > primal opt obj"
#     feasibility = ((log.dual_var[:σ]') .- log.dual_var[:δ] .- (log.dual_var[:λ] * d')) .<= c .+ 1e-6
#     @assert all(feasibility) "Dual feasibility violated"
#     return [Hyperplane(a_x, a_t, a_0)]
# end

# update y_k first, then update lambda (without f_x version)
function generate_cuts(data::Data, oracle::DualDecomposition, x::Vector{Float64})
    param = oracle.oracle_param; log = oracle.oracle_log
    y_k =log.pri_var; λ_k = log.dual_var[:λ]; best_obj = copy(log.best_obj); halving_value = copy(param.halving_value); α = copy(param.step_size)
    d, u = data.problem.demands, data.problem.capacities; c = data.problem.costs .* d'

    iter, γ = 1, 0
    while α >= param.stepsize_bound
        # Update y, compute residual and LB
        y_k, _, _, _, _ = update_y(y_k, λ_k, c, d, x)
        residual = [d'*y_k[i,:] - u[i]*x[i] for i in eachindex(u)]
        z_lb = sum(c .* y_k) + λ_k' * residual

        # Update step size
        α = halving_value * abs(param.obj_limit - z_lb)/norm(residual)^2 # Polyak's stepsize
        # α = 0.1/sqrt(iter) # diminishing

        # Update dual value
        λ_k = update_λ(λ_k, y_k, u, d, x, α)

        # Keep the best solution and obj
        if best_obj < z_lb 
            best_obj = z_lb
        else
            γ += 1
            if γ > param.γ_max
                halving_value = halving_value * 0.5
                γ = 0
            end
        end
        push!(log.α_set, α); push!(log.obj_set, z_lb)
        # push!(log.λ_k_norm_set, norm(λ_k)); push!(log.λ_k_diff_set, norm(opt_λ .+ λ_k)); push!(log.residual_norm_set, norm(residual)); push!(log.approx_obj_set, sum(c .* y_k))

        # Handle divergence of Polyak caused by not diminishing stepsize (do not immediately terminate)
        (best_obj > param.obj_limit && iter > param.max_iter) && (@info "best f: $(best_obj) > overestimated f: $(param.obj_limit)"; log.flag_exceed = true; break)

        iter += 1
    end
    _, sorted_indices, c_sorted, _, critical_facility_indices = update_y(y_k, λ_k, c, d, x)
    retrieve_dual_values(log, d, sorted_indices, c_sorted, critical_facility_indices)

    a_x = -log.dual_var[:λ].*u .- [sum(log.dual_var[:δ][i,:]) for i in eachindex(u)]
    a_t = [-1.0]
    a_0 = sum(log.dual_var[:σ]) 

    feasibility = ((log.dual_var[:σ]') .- log.dual_var[:δ] .- (log.dual_var[:λ] * d')) .<= c .+ 1e-6
    @assert all(feasibility) "Dual feasibility violated"
    return [Hyperplane(a_x, a_t, a_0)]
end

# Update lambda firt, then update y_k
# function generate_cuts(data::Data, oracle::DualDecomposition, x::Vector{Float64}, f_x; time_limit = 3600)
#     param = oracle.oracle_param; log = oracle.oracle_log
#     y_k =log.pri_var; λ_k = log.dual_var[:λ]; best_obj = copy(log.best_obj); halving_value = copy(param.halving_value); α = copy(param.step_size)
#     d, u = data.problem.demands, data.problem.capacities; c = data.problem.costs .* d'

#     iter, γ = 1, 0
#     sorted_indices = nothing; c_sorted = nothing; critical_facility_indices =nothing
#     @info param.obj_limit

#     while α >= param.stepsize_bound
#         # Update dual value
#         λ_k = update_λ(λ_k, y_k, u, d, x, α)

#         # Update y, compute residual and LB
#         y_k, sorted_indices, c_sorted, _, critical_facility_indices = update_y(y_k, λ_k, c, d, x)
#         residual = [d'*y_k[i,:] - u[i]*x[i] for i in eachindex(u)]
#         z_lb = sum(c .* y_k) + λ_k' * residual

#         # Update step size
#         α = halving_value * abs(param.obj_limit - z_lb)/norm(residual)^2 # Polyak's stepsize
#         # α = α/sqrt(iter) # diminishing

#         # Keep the best solution and obj
#         if best_obj < z_lb 
#             best_obj = z_lb
#         else
#             γ += 1
#             if γ > param.γ_max
#                 halving_value = halving_value * 0.5
#                 γ = 0
#             end
#         end
#         push!(log.α_set, α); push!(log.obj_set, z_lb)
#         # push!(log.λ_k_norm_set, norm(λ_k)); push!(log.λ_k_diff_set, norm(opt_λ .+ λ_k)); push!(log.residual_norm_set, norm(residual)); push!(log.approx_obj_set, sum(c .* y_k))

#         # Handle divergence of Polyak caused by not diminishing stepsize (do not immediately terminate)
#         (best_obj > param.obj_limit && iter > param.max_iter) && (@info "best f: $(best_obj) > overestimated f: $(param.obj_limit)"; log.flag_exceed = true; break)

#         iter += 1
#     end
#     retrieve_dual_values(log, d, sorted_indices, c_sorted, critical_facility_indices)

#     a_x = -log.dual_var[:λ].*u .- [sum(log.dual_var[:δ][i,:]) for i in eachindex(u)]
#     a_t = [-1.0]
#     a_0 = sum(log.dual_var[:σ]) 

#     @assert (sum(log.dual_var[:σ]) - dot(x, vec(sum(log.dual_var[:δ],dims = 2))) - dot(log.dual_var[:λ], u.*x) <= f_x +1e-3) "Dual obj > primal opt obj"
#     feasibility = ((log.dual_var[:σ]') .- log.dual_var[:δ] .- (log.dual_var[:λ] * d')) .<= c .+ 1e-6
#     @assert all(feasibility) "Dual feasibility violated"
#     return [Hyperplane(a_x, a_t, a_0)]
# end