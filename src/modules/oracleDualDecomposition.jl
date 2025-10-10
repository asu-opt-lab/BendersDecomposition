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
    vogel_time::Float64

    function DualDecompositionLog(data::Data; obj_set = [], α_set = [], λ_k_diff_set = [], residual_norm_set = [], num_solver_used = 0, prev_is_in_L = false, y_diff_set = [], fin_LB_set = [], vogel_time=Inf)
        m = data.dim_x; n = data.problem.n_customers
        pri_var =  zeros(Float64, m, n); dual_var = Dict(:λ => zeros(m), :σ => zeros(n), :δ => zeros(Float64, m, n))
        new(pri_var, dual_var, obj_set, α_set, λ_k_diff_set, residual_norm_set, num_solver_used, prev_is_in_L, y_diff_set, fin_LB_set, vogel_time)
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

# update order: y → λ 
function generate_cuts(oracle::DualDecomposition, x::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    data = oracle.data; param = oracle.oracle_param; log = oracle.oracle_log; typical_oracle = oracle.typical_oracle

    if log.prev_is_in_L || (length(log.fin_LB_set) >= param.stale_lim + 1 && all(isapprox.(log.fin_LB_set[end], log.fin_LB_set[end-param.stale_lim:end-1])))
        println("Exact solver")
        log.num_solver_used += 1; log.prev_is_in_L = false
        is_in_L, hyperplanes, f_x = generate_cuts(typical_oracle, x, t_value; time_limit = time_limit)
        return is_in_L, hyperplanes, f_x
    else
        println("Dual Decomposition")
        log.vogel_time = @elapsed param.obj_limit = vogel(data, x);

        # Just for comparison with optimal info
        _, _, f_x = generate_cuts(typical_oracle, x, t_value)
        opt_y = value.(typical_oracle.model[:y]); opt_λ = dual.(typical_oracle.model[:capacity]); @debug opt_λ

        # call necessary parameters
        y_k =log.pri_var; λ_k = log.dual_var[:λ]; halving_value = copy(param.halving_value); α = copy(param.step_size)
        d = data.problem.demands; u= data.problem.capacities; c = data.problem.costs .* d'

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

            # Update halving parameter
            γ += 1
            γ > param.γ_max && (halving_value = halving_value * 0.5; γ = 0)

            # Store necessary results
            push!(log.α_set, α); push!(log.obj_set, z_lb); # push!(log.residual_norm_set, norm(residual)); push!(log.λ_k_diff_set, norm(opt_λ .+ λ_k))

            iter += 1
        end
        y_k, sorted_indices, c_sorted, _, critical_facility_indices = update_y(y_k, λ_k, c, d, x)
        retrieve_dual_values(log, d, sorted_indices, c_sorted, critical_facility_indices)
        
        feasibility = ((log.dual_var[:σ]') .- log.dual_var[:δ] .- (log.dual_var[:λ] * d')) .<= c .+ 1e-6
        @assert all(feasibility) "Dual feasibility violated"

        dual_obj = sum(log.dual_var[:σ]) - dot(x, vec(sum(log.dual_var[:δ],dims = 2))) - dot(log.dual_var[:λ], u.*x)
        push!(log.fin_LB_set, dual_obj)

        log.prev_is_in_L = false
        if dual_obj >= t_value[1]/tol_normalize * (1 + 1e-6) + 1e-6/tol_normalize
            a_x = -log.dual_var[:λ].*u .- [sum(log.dual_var[:δ][i,:]) for i in eachindex(u)]
            a_t = [-1.0]
            a_0 = sum(log.dual_var[:σ]) 
            # return false, [Hyperplane(a_x, a_t, a_0)], [NaN]
            return false, [Hyperplane(a_x, a_t, a_0)], f_x
        else
            if !oracle.flag_bnb
                log.prev_is_in_L = true
                # return false, [Hyperplane(zeros(length(x)), [0.0], 0.0)], [NaN]
                return false, [Hyperplane(zeros(length(x)), [0.0], 0.0)], f_x
            else
                is_in_L, hyperplanes, f_x = generate_cuts(typical_oracle, x, t_value; time_limit = time_limit)
                println("lambda_DD:$(log.dual_var[:λ])")
                return is_in_L, hyperplanes, f_x
            end
        end
    end
end