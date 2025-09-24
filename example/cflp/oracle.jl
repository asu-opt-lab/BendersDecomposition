struct FacilityKnapsackInfo
    costs::Matrix{Float64}
    demands::Vector{Float64}
    capacity::Vector{Float64}
end

mutable struct CFLKnapsackOracleParam <: AbstractOracleParam
    atol::Float64
    rtol::Float64

    function CFLKnapsackOracleParam(; atol = 1e-9, rtol = 1e-9)
        new(atol, rtol)
    end
end

mutable struct CFLKnapsackOracle <: AbstractTypicalOracle
    oracle_param::CFLKnapsackOracleParam

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    facility_knapsack_info::FacilityKnapsackInfo

    function CFLKnapsackOracle(data::Data; 
                               scen_idx=-1, 
                               solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9),
                               oracle_param::CFLKnapsackOracleParam = CFLKnapsackOracleParam())
        @debug "Building classical oracle"
        model = Model()

        # Define coupling variables and constraints
        @variable(model, 0 <= x[1:data.dim_x] <= 1)
        @constraint(model, fix_x, x .== 0)

        facility_knapsack_info = scen_idx == -1 ? FacilityKnapsackInfo(data.problem.costs, data.problem.demands, data.problem.capacities) : FacilityKnapsackInfo(data.problem.costs, data.problem.demands[scen_idx], data.problem.capacities)

        assign_attributes!(model, solver_param)
        
        new(oracle_param, model, fix_x, facility_knapsack_info)
    end
    
    CFLKnapsackOracle() = new()
end

function generate_cuts(oracle::CFLKnapsackOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    optimize!(oracle.model)
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    status = dual_status(oracle.model)

    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)

        μ = dual.(oracle.model[:demand])
        a_t = [-1.0] 
        
        # Get facility knapsack info
        costs = oracle.facility_knapsack_info.costs
        demands = oracle.facility_knapsack_info.demands
        capacity = oracle.facility_knapsack_info.capacity

        # Calculate KP values for each facility
        KP_values = Vector{Float64}(undef, length(capacity))
        for i in 1:length(capacity)
            KP_values[i] = calculate_KP_value(costs[i,:], demands, capacity[i], μ)
        end

        a_x = KP_values # Vector{Float64}
        a_0 = sum(μ) 
        if sub_obj_val >= t_value[1] * (1 + oracle.oracle_param.rtol / tol_normalize)
            return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], deepcopy(t_value)
        end
        
    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(oracle.model)
            dual_sub_obj_val = dual_objective_value(oracle.model)
            @info "dual_sub_obj_val = $dual_sub_obj_val"
            
            a_x = dual.(oracle.fixed_x_constraints)
            a_t = [0.0]
            a_0 = dual_sub_obj_val - a_x' * x_value 
            if dual_sub_obj_val >= oracle.oracle_param.atol / tol_normalize
                return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
            else
                return true, [Hyperplane(a_x, a_t, a_0)], [Inf]
            end
        end
        
    else
        throw(UnexpectedModelStatusException("ClassicalOracle: $(status)"))
    end
end

function calculate_KP_value(costs::Vector{Float64}, demands::Vector{Float64}, capacity::Float64, μ::Vector{Float64})
    n = length(demands)
    
    # ratios = Vector{Tuple{Int,Float64}}(undef, n)
    ratios = [(i, (costs[i] * demands[i] - μ[i]) / demands[i]) for i in 1:n if (costs[i] * demands[i] - μ[i]) < 0]
    
    sort!(ratios, by=x->x[2])
    
    kp_value = 0.0
    remaining_capacity = capacity
    z = zeros(n)

    for (i, _) in ratios
        if remaining_capacity >= demands[i]
            kp_value += costs[i] * demands[i] - μ[i]
            remaining_capacity -= demands[i]
            z[i] = 1.0
        else
            fraction = remaining_capacity / demands[i]
            kp_value += (costs[i]*demands[i] - μ[i]) * fraction
            z[i] = fraction
            break
        end
    end

    return kp_value
end

# function northwest_corner(data, facilities::Vector{Float64})
#     assignments = zeros(Float64, data.dim_x, data.problem.n_customers)
#     opened = findall(x -> x == 1.0, facilities)
#     @assert !isempty(opened) "No facilities openend"

#     demands, capacities = copy(data.problem.demands), copy(data.problem.capacities)
#     @assert capacities' * facilities >= sum(demands) "Feasibility condition is violated, capa: $(capacities' * facilities), demand: $(sum(demands))"

#     fulfillment = zeros(length(demands))
#     i, j = 1, 1

#     # Do Northwest corner method
#     while i <= length(capacities) && j <= length(demands)
#         # Check if i ∉ opened facilities.
#         !in(i, opened) && (i +=1; continue)

#         # Assign demands to facilities               
#         rem = 1 - fulfillment[j]
#         take = min(rem, capacities[i] / demands[j])
#         assignments[i, j] = take
#         fulfillment[j] += take     
#         isapprox(fulfillment[j], 1.0, atol=1e-4) && (fulfillment[j] = 1.0)
#         capacities[i] -= take * demands[j]
#         isapprox(capacities[i], 0.0, atol=1e-4) && (capacities[i] = 0.0)

#         # Index update rule
#         capacities[i] == 0.0 && (i += 1)
#         fulfillment[j] == 1.0 && (j += 1)
#     end
#     return sum((data.problem.costs .* data.problem.demands') .* assignments)
# end

# mine original
function vogel(data, facilities::Vector{Float64})
    demands, capacities = copy(data.problem.demands), copy(data.problem.capacities)
    costs =  data.problem.costs .* demands'
    not_opened = findall(x -> x == 0.0, facilities)
    costs[not_opened,:] .= +Inf

    assignments = zeros(Float64, data.dim_x, data.problem.n_customers)
    row_diff = zeros(Float64, data.dim_x); col_diff = zeros(Float64, data.problem.n_customers)
    fulfillment = zeros(length(demands))

    while any(fulfillment .< 1 - 1e-6)
        for idx in eachindex(capacities)
            v = costs[idx,:]
            if all(isinf, v)
                row_diff[idx] = -Inf
            else
                small1, small2 = partialsort(v, 1:2)
                row_diff[idx] = small2 - small1
            end
        end

        for idx in eachindex(demands)
            v = costs[:, idx]
            if all(isinf, v)
                col_diff[idx] = -Inf
            else
                small1, small2 = partialsort(v, 1:2)
                col_diff[idx] = small2 - small1
            end
        end

        i = nothing; j = nothing
        if maximum(row_diff) >= maximum(col_diff)
            i = argmax(row_diff)
            j = argmin(costs[i,:])
        else
            j = argmax(col_diff)
            i = argmin(costs[:,j])
        end

        rem = 1 - fulfillment[j]
        take = min(rem, capacities[i] / demands[j])
        assignments[i, j] = take
        fulfillment[j] += take     
        capacities[i] -= take * demands[j]

        isapprox(capacities[i], 0.0; atol=1e-6) && (costs[i,:] .= Inf)
        isapprox(fulfillment[j], 1.0; atol=1e-6) && (costs[:, j] .= Inf)
    end
    return sum((data.problem.costs.* demands') .* assignments)
end

# mine refined, but takes much time than original.
# function vogel(data, facilities::Vector{Float64})
#     demands, capacities = copy(data.problem.demands), copy(data.problem.capacities)
#     costs =  data.problem.costs .* demands'
#     not_opened = findall(x -> x == 0.0, facilities)
#     costs_not_opened = @view costs[not_opened,:]; costs_not_opened .= +Inf

#     assignments = zeros(Float64, data.dim_x, data.problem.n_customers)
#     row_diff = zeros(Float64, data.dim_x); col_diff = zeros(Float64, data.problem.n_customers)
#     fulfillment = zeros(length(demands))

#     while any(fulfillment .< 1 - 1e-6)
#         @info sum(fulfillment)
#         for idx in eachindex(capacities)
#             costs_idx = @view costs[idx,:]
#             small = Inf; small2 = Inf; anyfinite = false
#             for j in eachindex(demands)
#                 cost = costs_idx[j]
#                 if isfinite(cost)
#                     anyfinite = true
#                     if cost < small
#                         small2 = small
#                         small = cost
#                     elseif cost < small2
#                         small2 = cost
#                     end
#                 end
#             end
#             # If all elements = +Inf, set row_diff = -Inf so that the row is not selected. Else, we compute the difference.
#             row_diff[idx] = anyfinite ? (isfinite(small2) ? (small2 - small) : Inf) : -Inf
#         end

#         for idx in eachindex(demands)
#             costs_idx = @view costs[:, idx]
#             small = Inf; small2 = Inf; anyfinite = false
#             for i in eachindex(capacities)
#                 cost = costs_idx[i]
#                 if isfinite(cost)
#                     anyfinite = true
#                     if cost < small
#                         small2 = small
#                         small = cost
#                     elseif cost < small2
#                         small2 = cost
#                     end
#                 end
#             end
#             # If all elements = +Inf, set col_diff = -Inf so that the row is not selected. Else, we compute the difference.
#             col_diff[idx] = anyfinite ? (isfinite(small2) ? (small2 - small) : Inf) : -Inf
#         end

#         i = nothing; j = nothing
#         if maximum(row_diff) >= maximum(col_diff)
#             i = argmax(row_diff); j = argmin(@view costs[i,:])
#         else
#             j = argmax(col_diff); i = argmin(@view costs[:,j])
#         end

#         rem = 1 - fulfillment[j]
#         take = min(rem, capacities[i] / demands[j])
#         assignments[i, j] = take; fulfillment[j] += take; capacities[i] -= take * demands[j]

#         costs_i = @view costs[i,:]; costs_j = @view costs[:, j]
#         isapprox(capacities[i], 0.0; atol=1e-6) && (costs_i .= Inf)
#         isapprox(fulfillment[j], 1.0; atol=1e-6) && (costs_j .= Inf)
#     end
#     return sum((data.problem.costs.* demands') .* assignments)
# end

# GPT understanding
# function vogel(data, facilities::Vector{Float64})
#     # --- problem data ---
#     demands = copy(data.problem.demands); capacities = copy(data.problem.capacities)
#     open_mask = facilities .> 0.0

#     # working costs for selection (mask with Inf); assignments store FRACTIONS of demand
#     costs = data.problem.costs .* demands'           
#     costs[.!open_mask, :] .= Inf
#     m, n = size(costs)

#     assignments = zeros(Float64, m, n)
#     fulfillment = zeros(Float64, n)

#     # row/cols should be considered
#     row_active = copy(open_mask)
#     col_active = trues(n)

#     # penalties and argmins
#     row_diff   = fill(-Inf, m);  row_argmin = fill(0, m) # cheapest customer index of each row
#     col_diff   = fill(-Inf, n);  col_argmin = fill(0, n) # cheapest facility index of each column

#     # reverse maps: who currently points to whom (for local fixes)
#     rows_by_best_col = [Int[] for _ in 1:n] # what are the facility indices that consider j as the best? -> if j is satisfied, corresponding facilities that consider j as the best are influenced. the cost should be recomputed.
#     cols_by_best_row = [Int[] for _ in 1:m] # what are the customer indices that consider i as the best?

#     # ---- helpers: recompute ONE row/column, touching only active counterpart indices ----
#     function recompute_row!(i)
#         prev = row_argmin[i] # previous cheapest customer for the facility
#         best, best2, jbest = Inf, Inf, 0
#         @inbounds for j in 1:n
#             col_active[j] || continue # check jth col should be considered
#             c = costs[i, j]; isfinite(c) || continue # +Inf value should not be considered
#             if c < best
#                 best2 = best; best = c; jbest = j
#             elseif c < best2
#                 best2 = c
#             end
#         end
#         row_argmin[i] = jbest # current cheapest customer
#         # If all elements are Inf, we cannot use that row
#         row_diff[i]   = (isfinite(best) && isfinite(best2)) ? (best2 - best) : (isfinite(best) ? Inf : -Inf)
#         if prev != 0
#             v = rows_by_best_col[prev] # facilities that considered jth customer as the cheapest
#             for k in eachindex(v); v[k] == i && (deleteat!(v, k); break); end
#         end
#         jbest != 0 && push!(rows_by_best_col[jbest], i)
#         return nothing
#     end

#     function recompute_col!(j)
#         prev = col_argmin[j]
#         best, best2, ibest = Inf, Inf, 0
#         @inbounds for i in 1:m
#             row_active[i] || continue
#             c = costs[i, j]; isfinite(c) || continue
#             if c < best
#                 best2 = best; best = c; ibest = i
#             elseif c < best2
#                 best2 = c
#             end
#         end
#         col_argmin[j] = ibest
#         col_diff[j]   = (isfinite(best) && isfinite(best2)) ? (best2 - best) :
#                         (isfinite(best) ? Inf : -Inf)
#         if prev != 0
#             v = cols_by_best_row[prev]
#             for k in eachindex(v); v[k] == j && (deleteat!(v, k); break); end
#         end
#         ibest != 0 && push!(cols_by_best_row[ibest], j)
#         return nothing
#     end

#     # ---- initial penalties (one pass over matrix) ----
#     for i in 1:m
#         row_active[i] && recompute_row!(i)
#     end
#     for j in 1:n
#         col_active[j] && recompute_col!(j)
#     end

#     tol = 1e-6
#     while any(fulfillment .< 1 - tol)
#         # --- pick a valid line and cell (refresh lazily if stale) ---
#         i, j = 0, 0
#         while true
#             i_star = argmax(row_diff)
#             j_star = argmax(col_diff)
#             r_pen  = row_diff[i_star]
#             c_pen  = col_diff[j_star]

#             if r_pen == -Inf && c_pen == -Inf
#                 error("No feasible cell remains, yet some demand is unmet.")
#             end

#             # need to assign customer to the facility first when r_pen >= c_pen
#             if r_pen >= c_pen
#                 # check ① facility closed, ② column is deactivated, ③ cost is Inf -> update the information because current information is unusable.
#                 if row_argmin[i_star] == 0 || !col_active[row_argmin[i_star]] || !isfinite(costs[i_star, row_argmin[i_star]])
#                     row_active[i_star] && recompute_row!(i_star)
#                     row_active[i_star] || (row_diff[i_star] = -Inf)
#                     continue
#                 end
#                 i = i_star
#                 j = row_argmin[i]
#                 break
#             else
#                 # ensure column is valid
#                 if col_argmin[j_star] == 0 || !row_active[col_argmin[j_star]] || !isfinite(costs[col_argmin[j_star], j_star])
#                     col_active[j_star] && recompute_col!(j_star)
#                     col_active[j_star] || (col_diff[j_star] = -Inf)
#                     continue
#                 end
#                 j = j_star
#                 i = col_argmin[j]
#                 break
#             end
#         end

#         # allocate (fractions)
#         rem  = 1 - fulfillment[j]
#         take = min(rem, capacities[i] / demands[j])
#         assignments[i, j] = take
#         fulfillment[j]   += take
#         capacities[i]    -= take * demands[j]

#         # deactivate and locally fix penalties only where needed
#         if capacities[i] <= 1e-6
#             row_active[i] = false
#             costs[i, :]   .= Inf
#             for jj in cols_by_best_row[i]
#                 col_active[jj] && recompute_col!(jj)
#             end
#             empty!(cols_by_best_row[i])
#             row_diff[i] = -Inf
#         end
#         if fulfillment[j] >= 1 - tol
#             col_active[j] = false
#             costs[:, j]   .= Inf
#             for ii in rows_by_best_col[j]
#                 row_active[ii] && recompute_row!(ii)
#             end
#             empty!(rows_by_best_col[j])
#             col_diff[j] = -Inf
#         end
#     end
#     return sum((data.problem.costs .* demands') .* assignments)
# end

