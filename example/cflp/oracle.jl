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

function northwest_corner(data, facilities::Vector{Float64})
    assignments = zeros(Float64, data.dim_x, data.problem.n_customers)
    opened = findall(x -> x == 1.0, facilities)
    @assert !isempty(opened) "No facilities openend"

    demands, capacities = copy(data.problem.demands), copy(data.problem.capacities)
    @assert capacities' * facilities >= sum(demands) "Feasibility condition is violated, capa: $(capacities' * facilities), demand: $(sum(demands))"

    fulfillment = zeros(length(demands))
    i, j = 1, 1

    # Do Northwest corner method
    while i <= length(capacities) && j <= length(demands)
        # Check if i ∉ opened facilities.
        !in(i, opened) && (i +=1; continue)

        # Assign demands to facilities               
        rem = 1 - fulfillment[j]
        take = min(rem, capacities[i] / demands[j])
        assignments[i, j] = take
        fulfillment[j] += take     
        isapprox(fulfillment[j], 1.0, atol=1e-4) && (fulfillment[j] = 1.0)
        capacities[i] -= take * demands[j]
        isapprox(capacities[i], 0.0, atol=1e-4) && (capacities[i] = 0.0)

        # Index update rule
        capacities[i] == 0.0 && (i += 1)
        fulfillment[j] == 1.0 && (j += 1)
    end
    return sum((data.problem.costs .* data.problem.demands') .* assignments)
end

# mine ver.1 (update full)
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
    @assert isapprox(sum(fulfillment), data.problem.n_customers; atol=1e-6) "demand not fulfilled: $(sum(fulfillment)), num customers: $(data.problem.n_customers))"
    @assert all(capacities .>= -1e-6) "capacity exceeds: $(capacities)"
    return sum((data.problem.costs.* demands') .* assignments)
end

# mine ver.2 (lazy update) 
# function vogel(data, facilities::Vector{Float64})
#     demands = copy(data.problem.demands); capacities = copy(data.problem.capacities); costs =  data.problem.costs .* demands'
#     closed = findall(x -> x == 0.0, facilities); costs[closed,:] .= Inf
#     m,n = size(costs)

#     assignments = zeros(Float64, m, n); fulfillment = zeros(Float64, n) 
    
#     row_regret = fill(-Inf, m); col_regret = fill(-Inf, n) # cheapest - second cheapest
#     map_best_rtc = [Int[] for _ in 1:n]; map_best_ctr = [Int[] for _ in 1:m] # rows/cols regard col/row as the cheapest 
#     map_best2_rtc = [Int[] for _ in 1:n]; map_best2_ctr = [Int[] for _ in 1:m] # rows/cols regard col/row as the 2nd cheapest

#     curr_best_col_of_row = fill(0, m); curr_best_row_of_col = fill(0, n) # Cache current cheapest column/row of the row/column
#     curr_best2_col_of_row = fill(0, m); curr_best2_row_of_col = fill(0, n) # Cache current 2nd cheapest column/row of the row/column

#     costs_idx_sorted_rows = [sortperm(@view costs[i, :]) for i in 1:m]; costs_idx_sorted_cols = [sortperm(@view costs[:, j]) for j in 1:n]

#     function recompute_rows!(i, flag)
#         @inbounds begin
#             if flag === true && !isempty(costs_idx_sorted_rows[i])
#                 popfirst!(costs_idx_sorted_rows[i])
#             elseif flag === false && length(costs_idx_sorted_rows[i]) >= 2
#                 deleteat!(costs_idx_sorted_rows[i], 2)
#             end

#             while !isempty(costs_idx_sorted_rows[i]) && !isfinite(costs[i, costs_idx_sorted_rows[i][1]])
#                 popfirst!(costs_idx_sorted_rows[i])
#             end
#             if length(costs_idx_sorted_rows[i]) >= 2 && !isfinite(costs[i, costs_idx_sorted_rows[i][2]])
#                 k = 3
#                 while k <= length(costs_idx_sorted_rows[i]) && !isfinite(costs[i, costs_idx_sorted_rows[i][k]])
#                     k += 1
#                 end
#                 if k <= length(costs_idx_sorted_rows[i])
#                     idx = costs_idx_sorted_rows[i][k]
#                     deleteat!(costs_idx_sorted_rows[i], k)
#                     insert!(costs_idx_sorted_rows[i], 2, idx)
#                 else
#                     costs_idx_sorted_rows[i] = costs_idx_sorted_rows[i][1:1]
#                 end
#             end

#             if curr_best_col_of_row[i] > 0
#                 c = curr_best_col_of_row[i]
#                 filter!(x -> x != i, map_best_rtc[c])
#                 curr_best_col_of_row[i] = 0
#             end
#             if curr_best2_col_of_row[i] > 0
#                 c2 = curr_best2_col_of_row[i]
#                 filter!(x -> x != i, map_best2_rtc[c2])
#                 curr_best2_col_of_row[i] = 0
#             end

#             if isempty(costs_idx_sorted_rows[i])
#                 row_regret[i] = -Inf
#                 return
#             end

#             best_j = costs_idx_sorted_rows[i][1]
#             push!(map_best_rtc[best_j], i)
#             curr_best_col_of_row[i] = best_j

#             if length(costs_idx_sorted_rows[i]) > 1
#                 best2_j = costs_idx_sorted_rows[i][2]
#                 push!(map_best2_rtc[best2_j], i)
#                 curr_best2_col_of_row[i] = best2_j
#                 row_regret[i] = costs[i, best2_j] - costs[i, best_j]
#             else
#                 row_regret[i] = Inf
#             end
#         end
#     end

#     function recompute_cols!(j, flag)
#         @inbounds begin
#             if flag === true && !isempty(costs_idx_sorted_cols[j])
#                 popfirst!(costs_idx_sorted_cols[j])
#             elseif flag === false && length(costs_idx_sorted_cols[j]) >= 2
#                 deleteat!(costs_idx_sorted_cols[j], 2)
#             end

#             while !isempty(costs_idx_sorted_cols[j]) && !isfinite(costs[costs_idx_sorted_cols[j][1], j])
#                 popfirst!(costs_idx_sorted_cols[j])
#             end
#             if length(costs_idx_sorted_cols[j]) >= 2 && !isfinite(costs[costs_idx_sorted_cols[j][2], j])
#                 k = 3
#                 while k <= length(costs_idx_sorted_cols[j]) && !isfinite(costs[costs_idx_sorted_cols[j][k], j])
#                     k += 1
#                 end
#                 if k <= length(costs_idx_sorted_cols[j])
#                     idx = costs_idx_sorted_cols[j][k]
#                     deleteat!(costs_idx_sorted_cols[j], k)
#                     insert!(costs_idx_sorted_cols[j], 2, idx)
#                 else
#                     costs_idx_sorted_cols[j] = costs_idx_sorted_cols[j][1:1]
#                 end
#             end

#             if curr_best_row_of_col[j] > 0
#                 r = curr_best_row_of_col[j]
#                 filter!(x -> x != j, map_best_ctr[r])
#                 curr_best_row_of_col[j] = 0
#             end
#             if curr_best2_row_of_col[j] > 0
#                 r2 = curr_best2_row_of_col[j]
#                 filter!(x -> x != j, map_best2_ctr[r2])
#                 curr_best2_row_of_col[j] = 0
#             end

#             if isempty(costs_idx_sorted_cols[j])
#                 col_regret[j] = -Inf
#                 return
#             end

#             best_i = costs_idx_sorted_cols[j][1]
#             push!(map_best_ctr[best_i], j)
#             curr_best_row_of_col[j] = best_i

#             if length(costs_idx_sorted_cols[j]) > 1
#                 best2_i = costs_idx_sorted_cols[j][2]
#                 push!(map_best2_ctr[best2_i], j)
#                 curr_best2_row_of_col[j] = best2_i
#                 col_regret[j] = costs[best2_i, j] - costs[best_i, j]
#             else
#                 col_regret[j] = Inf
#             end
#         end
#     end
#     # Compute row_regret
#     @inbounds for i in 1:m
#         (i in closed) && continue
#         recompute_rows!(i, nothing)
#     end

#     # Compute col_regret
#     @inbounds for j in 1:n
#         (demands[j] == 0) && continue
#         recompute_cols!(j, nothing)
#     end

#     i = 0; j = 0
#     while any(fulfillment .< 1 - 1e-6)
#         @inbounds begin
#             i = argmax(row_regret); r_pen = row_regret[i]
#             j = argmax(col_regret); c_pen = col_regret[j]

#             if r_pen > c_pen
#                 if isempty(costs_idx_sorted_rows[i])
#                     row_regret[i] = -Inf; continue
#                 end
#                 j = costs_idx_sorted_rows[i][1]
#                 if !isfinite(costs[i, j])
#                     recompute_rows!(i, true); continue
#                 end
#             else
#                 if isempty(costs_idx_sorted_cols[j])
#                     col_regret[j] = -Inf; continue
#                 end
#                 i = costs_idx_sorted_cols[j][1]
#                 if !isfinite(costs[i, j])
#                     recompute_cols!(j, true); continue
#                 end
#             end

#             rem = 1 - fulfillment[j]
#             take = min(rem, capacities[i] / demands[j])
#             assignments[i, j] = take; fulfillment[j] += take; capacities[i] -= take * demands[j]

#             if capacities[i] <= 1e-6
#                 costs[i, :] .= Inf
#                 for col in map_best_ctr[i]
#                     recompute_cols!(col, true)
#                 end
#                 empty!(map_best_ctr[i])
#                 for col in map_best2_ctr[i]
#                     recompute_cols!(col, false)
#                 end
#                 empty!(map_best2_ctr[i])
#                 row_regret[i] = -Inf
#             end
            
#             if fulfillment[j] >= 1 - 1e-6
#                 costs[:, j] .= Inf
#                 for row in map_best_rtc[j]
#                     recompute_rows!(row, true)
#                 end
#                 empty!(map_best_rtc[j])
#                 for row in map_best2_rtc[j]
#                     recompute_rows!(row, false)
#                 end
#                 empty!(map_best2_rtc[j])
#                 col_regret[j] = -Inf
#             end
#         end
#     end
#     @assert isapprox(sum(fulfillment), n; atol=1e-6) "demand not fulfilled: $(sum(fulfillment)), num customers: $n"
#     @assert all(capacities .>= -1e-6) "capacity exceeds: $(capacities)"
#     return sum((data.problem.costs.* demands') .* assignments)
# end

# vogel fastest
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
#     @assert isapprox(sum(fulfillment), n; atol=1e-6) "demand not fulfilled: $(sum(fulfillment)), num customers: $n"
#     @assert all(capacities .>= -1e-6) "capacity exceeds: $(capacities)"
#     return sum((data.problem.costs .* demands') .* assignments)
# end