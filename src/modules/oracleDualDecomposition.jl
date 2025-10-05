export DualDecompositionParam, DualDecompositionLog, DualDecomposition

include("$(dirname(dirname(@__DIR__)))/example/cflp/oracle.jl")

mutable struct DualDecompositionParam <: AbstractOracleParam
    stale_lim::Int64
    max_iter::Float64
    obj_limit::Float64
    step_size::Float64
    halving_value::Float64
    γ_max::Float64
    stepsize_bound::Float64

    function DualDecompositionParam(; stale_lim = 3, max_iter = 1000, obj_limit = +Inf, step_size = 0.1, halving_value = 2.0, γ_max = max_iter*0.1, stepsize_bound = 1e-4)
        new(stale_lim, max_iter, obj_limit, step_size, halving_value, γ_max, stepsize_bound)
    end
end

mutable struct DualDecompositionLog <: AbstractOracleLog
    pri_var::Matrix{Float64}
    dual_var::Dict{Symbol, Array{Float64}}
    obj_set::Vector{Any}
    α_set::Vector{Any}
    λ_k_diff_set::Vector{Any}
    residual_norm_set::Vector{Any}
    num_solver_used::Float64
    prev_is_in_L::Bool
    y_diff_set::Vector{Any}
    fin_LB_set::Vector{Any}

    function DualDecompositionLog(data::Data; obj_set = [], α_set = [], λ_k_diff_set = [], residual_norm_set = [], num_solver_used = 0, prev_is_in_L = false, y_diff_set = [], fin_LB_set = [])
        m = data.dim_x; n = data.problem.n_customers
        pri_var =  zeros(Float64, m, n); dual_var = Dict(:λ => zeros(m), :σ => zeros(n), :δ => zeros(Float64, m, n))
        new(pri_var, dual_var, obj_set, α_set, λ_k_diff_set, residual_norm_set, num_solver_used, prev_is_in_L, y_diff_set, fin_LB_set)
    end
end

mutable struct DualDecomposition <: AbstractTypicalOracle
    data::Data
    typical_oracle::AbstractOracle
    oracle_param::DualDecompositionParam
    oracle_log::DualDecompositionLog

    function DualDecomposition(data::Data; typical_oracle::AbstractOracle = ClassicalOracle(data),
        oracle_param::DualDecompositionParam = DualDecompositionParam(), oracle_log::DualDecompositionLog = DualDecompositionLog(data))
        @debug "Building Dual Decomposition Oracle"
        new(data, typical_oracle, oracle_param, oracle_log)
    end

    DualDecomposition() = new()
end

# update order: y → λ 
# function generate_cuts(oracle::DualDecomposition, x::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
#     data = oracle.data; param = oracle.oracle_param; log = oracle.oracle_log; typical_oracle = oracle.typical_oracle

#     if log.prev_is_in_L || (length(log.fin_LB_set) >= param.stale_lim + 1 && all(isapprox.(log.fin_LB_set[end], log.fin_LB_set[end-param.stale_lim:end-1])))
#         @debug "Exact solver"
#         log.num_solver_used += 1; log.prev_is_in_L = false
#         is_in_L, hyperplanes, f_x = generate_cuts(typical_oracle, x, t_value; time_limit = time_limit)
#         return is_in_L, hyperplanes, f_x
#     else
#         @debug "Dual Decomposition"
#         t = @elapsed param.obj_limit = vogel(data, x); @debug "Vogel obj: $(param.obj_limit), time_spent: $t"

#         # Just for comparison with optimal info
#         # _, _, f_x = generate_cuts(typical_oracle, x, t_value)
#         # opt_y = value.(typical_oracle.model[:y]); opt_λ = dual.(typical_oracle.model[:capacity]); @info opt_λ

#         # call necessary parameters
#         y_k =log.pri_var; λ_k = log.dual_var[:λ]; halving_value = copy(param.halving_value); α = copy(param.step_size)
#         d = data.problem.demands; u= data.problem.capacities; c = data.problem.costs .* d'

#         iter, γ = 1, 0
#         while α >= param.stepsize_bound
#             # Update y, compute residual and LB
#             y_k, _, _, _, _ = update_y(y_k, λ_k, c, d, x)
#             residual = [d'*y_k[i,:] - u[i]*x[i] for i in eachindex(u)]
#             z_lb = sum(c .* y_k) + λ_k' * residual

#             # Update step size
#             α = halving_value * abs(param.obj_limit - z_lb)/norm(residual)^2 # Polyak's stepsize
#             # α = 0.1/sqrt(iter) # diminishing

#             # Update dual value
#             λ_k = update_λ(λ_k, y_k, u, d, x, α)

#             # Update halving parameter
#             γ += 1
#             γ > param.γ_max && (halving_value = halving_value * 0.5; γ = 0)

#             # Store necessary results
#             push!(log.α_set, α); push!(log.obj_set, z_lb); # push!(log.residual_norm_set, norm(residual)); push!(log.λ_k_diff_set, norm(opt_λ .+ λ_k))

#             iter += 1
#         end
#         _, sorted_indices, c_sorted, _, critical_facility_indices = update_y(y_k, λ_k, c, d, x)
#         retrieve_dual_values(log, d, sorted_indices, c_sorted, critical_facility_indices)
        
#         feasibility = ((log.dual_var[:σ]') .- log.dual_var[:δ] .- (log.dual_var[:λ] * d')) .<= c .+ 1e-6
#         @assert all(feasibility) "Dual feasibility violated"


#         # gen_figure(log, length(log.fin_LB_set)); empty!(log.α_set); empty!(log.λ_k_diff_set) # Draw figure
#         @debug "λ^{approx}:$(log.dual_var[:λ])"

#         dual_obj = sum(log.dual_var[:σ]) - dot(x, vec(sum(log.dual_var[:δ],dims = 2))) - dot(log.dual_var[:λ], u.*x)
#         push!(log.fin_LB_set, dual_obj)
#         # @debug "primal: $(f_x[1]), dual: $dual_obj"
#         if dual_obj >= t_value[1]/tol_normalize * (1 + 1e-6) + 1e-6/tol_normalize
#             a_x = -log.dual_var[:λ].*u .- [sum(log.dual_var[:δ][i,:]) for i in eachindex(u)]
#             a_t = [-1.0]
#             a_0 = sum(log.dual_var[:σ]) 

#             log.prev_is_in_L = false
#             # log.pri_var = zeros(Float64, data.dim_x, data.problem.n_customers); log.dual_var[:λ] = zeros(data.dim_x) # reset y, λ
#             return false, [Hyperplane(a_x, a_t, a_0)], [NaN]
#         else
#             log.prev_is_in_L = typeof(typical_oracle) == EmptyOracle ? false : true # should not set true in BnB, will not use DD forever
#             # log.pri_var = zeros(Float64, data.dim_x, data.problem.n_customers); log.dual_var[:λ] = zeros(data.dim_x) # reset y, λ
#             return false, [Hyperplane(zeros(length(x)), [0.0], 0.0)], [NaN]
#         end
#     end
# end

# update order: λ → y
function generate_cuts(oracle::DualDecomposition, x::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    data = oracle.data; param = oracle.oracle_param; log = oracle.oracle_log; typical_oracle = oracle.typical_oracle

    if log.prev_is_in_L || (length(log.fin_LB_set) >= param.stale_lim + 1 && all(isapprox.(log.fin_LB_set[end], log.fin_LB_set[end-param.stale_lim:end-1])))
        @debug "Exact solver"
        log.num_solver_used += 1; log.prev_is_in_L = false
        is_in_L, hyperplanes, f_x = generate_cuts(typical_oracle, x, t_value; time_limit = time_limit)
        return is_in_L, hyperplanes, f_x
    else
        @debug "Dual Decomposition"
        t = @elapsed param.obj_limit = vogel(data, x); @debug (param.obj_limit, t)

        # Just for comparison with optimal info
        # _, _, f_x = generate_cuts(typical_oracle, x, t_value)
        # opt_y = value.(typical_oracle.model[:y]); opt_λ = dual.(typical_oracle.model[:capacity]); @debug opt_λ

        # call necessary parameters
        y_k =log.pri_var; λ_k = log.dual_var[:λ]; halving_value = copy(param.halving_value); α = copy(param.step_size)
        d = data.problem.demands; u= data.problem.capacities; c = data.problem.costs .* d'

        iter, γ = 1, 0
        sorted_indices = nothing; c_sorted = nothing; critical_facility_indices =nothing
        while α >= param.stepsize_bound
            # Update dual value
            λ_k = update_λ(λ_k, y_k, u, d, x, α)
            # y_i = @view y_k[7,:]
            # iter <= 5 && @info "first lambda: $(λ_k[7]), y_i_sum: $(sum(y_i)), α: $α, slack: $(dot(d, y_i) - u[7]*x[7])"

            # Update y, compute residual and LB
            y_k, sorted_indices, c_sorted, _, critical_facility_indices = update_y(y_k, λ_k, c, d, x)
            residual = [d'*y_k[i,:] - u[i]*x[i] for i in eachindex(u)]
            z_lb = sum(c .* y_k) + λ_k' * residual

            # Update step size
            α = halving_value * abs(param.obj_limit - z_lb)/norm(residual)^2 # Polyak's stepsize
            # α = 0.1/sqrt(iter) # diminishing

            # Update halving parameter
            γ += 1
            γ > param.γ_max && (halving_value = halving_value * 0.5; γ = 0)

            # Store necessary results
            # push!(log.α_set, α); push!(log.obj_set, z_lb); push!(log.residual_norm_set, norm(residual)); push!(log.λ_k_diff_set, norm(opt_λ .+ λ_k))

            iter += 1
        end
        retrieve_dual_values(log, d, sorted_indices, c_sorted, critical_facility_indices)

        feasibility = ((log.dual_var[:σ]') .- log.dual_var[:δ] .- (λ_k* d')) .<= c .+ 1e-6
        @assert all(feasibility) "Dual feasibility violated"

        # gen_figure(log, length(log.fin_LB_set)); empty!(log.α_set); empty!(log.λ_k_diff_set) # Draw figure
        @debug log.dual_var[:λ]

        dual_obj = sum(log.dual_var[:σ]) - dot(x, vec(sum(log.dual_var[:δ],dims = 2))) - dot(log.dual_var[:λ], u.*x)
        push!(log.fin_LB_set, dual_obj)
        # @debug "primal: $(f_x[1]), dual: $dual_obj"
        if dual_obj >= t_value[1]/tol_normalize * (1 + 1e-6) + 1e-6/tol_normalize
            a_x = -log.dual_var[:λ].*u .- [sum(log.dual_var[:δ][i,:]) for i in eachindex(u)]
            a_t = [-1.0]
            a_0 = sum(log.dual_var[:σ]) 

            log.prev_is_in_L = false
            log.dual_var[:λ] = zeros(data.dim_x); log.pri_var = zeros(Float64, data.dim_x, data.problem.n_customers) # doing this make invalid cut
            return false, [Hyperplane(a_x, a_t, a_0)], [NaN]
        else
            log.prev_is_in_L = typeof(typical_oracle) == EmptyOracle ? false : true # should not set true in BnB, will not use DD forever
            log.dual_var[:λ] = zeros(data.dim_x); log.pri_var = zeros(Float64, data.dim_x, data.problem.n_customers) # doing this make invalid cut
            return false, [Hyperplane(zeros(length(x)), [0.0], 0.0)], [NaN]
        end
    end
end