function generate_cuts(env::BendersEnv, ::KnapsackCut)
    (μ, KP_values, coeff_t), sub_obj_val = generate_cut_coefficients(env.sub, env.master.x_value, KnapsackCut())
    cut = @expression(env.master.model, 
        coeff_t * env.master.var[:t] + μ_term(env.sub, μ) + dot(KP_values, env.master.var[:x]))
    return cut, sub_obj_val
end

function generate_cuts(env::BendersEnv, ::KnapsackCut, scenario::Int)
    (μ, KP_values, coeff_t), sub_obj_val = generate_cut_coefficients(env.sub.sub_problems[scenario], env.master.x_value, KnapsackCut())
    cut = @expression(env.master.model, 
        coeff_t * env.master.var[:t][scenario] + μ_term(env.sub, μ) + dot(KP_values, env.master.var[:x]))

    return cut, sub_obj_val
end

μ_term(sub::AbstractSubProblem, μ::Float64) = μ
μ_term(sub::Union{KnapsackCFLPSubProblem, KnapsackSCFLPSubProblem}, μ::Vector{Float64}) = sum(μ)
μ_term(sub::KnapsackMCNDPSubProblem, μ::Matrix{Float64}) = sum(μ .* sub.b_iv)


function generate_cut_coefficients(sub::KnapsackCFLPSubProblem, x_value::Vector{Float64}, ::KnapsackCut)
    status = dual_status(sub.model)
    if status == FEASIBLE_POINT
        subObjVal = objective_value(sub.model)
        μ = dual.(sub.demand_constraints)

        # Get facility knapsack info
        costs = sub.facility_knapsack_info.costs
        demands = sub.facility_knapsack_info.demands
        capacity = sub.facility_knapsack_info.capacity

        # Calculate KP values for each facility
        KP_values = Vector{Float64}(undef, length(capacity))
        for i in 1:length(capacity)
            KP_values[i] = calculate_KP_value(costs[i,:], demands, capacity[i], μ)
        end

        return (μ, KP_values, -1.0), subObjVal

    elseif status == INFEASIBILITY_CERTIFICATE
        @debug "Infeasible subproblem"
        if has_duals(sub.model)
            coefficients_x = dual.(sub.fixed_x_constraints)
            constant_term = dot(dual.(sub.other_constraints), normalized_rhs.(sub.other_constraints))
            return (constant_term, coefficients_x, 0.0), Inf
        else
            throw(ErrorException("Infeasible subproblem has no dual solution"))
        end
        
    else
        @info "Unexpected dual status: $status"
        throw(ErrorException("Unexpected dual status"))
    end
end

function generate_cut_coefficients(sub::KnapsackMCNDPSubProblem, x_value::Vector{Float64}, ::KnapsackCut)
    status = dual_status(sub.model)
    if status == FEASIBLE_POINT
        subObjVal = objective_value(sub.model)
        μ = dual.(sub.demand_constraints)
        _μ = Matrix{Float64}(undef, sub.data.num_commodities, sub.data.num_arcs)
        for c in 1:sub.data.num_commodities
            for a in 1:sub.data.num_arcs
                from_node = sub.data.arcs[a][1]
                to_node = sub.data.arcs[a][2]
                _μ[c,a] = μ[c,from_node] - μ[c,to_node]
            end
        end

        costs = sub.data.variable_costs
        demands = sub.data.demands
        capacity = sub.data.capacities

        # Calculate KP values for each facility
        KP_values = Vector{Float64}(undef, length(capacity))
        for j in 1:length(capacity)
            KP_values[j] = calculate_KP_value(
                costs[j]*ones(sub.data.num_commodities), 
                [d[3] for d in demands], 
                capacity[j], 
                _μ[:,j])
        end

        return (μ, KP_values, -1.0), subObjVal

    elseif status == INFEASIBILITY_CERTIFICATE
        @debug "Infeasible subproblem"
        if has_duals(sub.model)
            μ = dual.(sub.demand_constraints)
            _μ = Matrix{Float64}(undef, sub.data.num_commodities, sub.data.num_arcs)
            for c in 1:sub.data.num_commodities
                for a in 1:sub.data.num_arcs
                    from_node = sub.data.arcs[a][1]
                    to_node = sub.data.arcs[a][2]
                    _μ[c,a] = μ[c,from_node] - μ[c,to_node]
                end
            end

            costs = sub.data.variable_costs
            demands = sub.data.demands
            capacity = sub.data.capacities

            # Calculate KP values for each facility
            KP_values = Vector{Float64}(undef, length(capacity))
            for j in 1:length(capacity)
                KP_values[j] = calculate_KP_value_infeasible(
                    costs[j]*ones(sub.data.num_commodities), 
                    [d[3] for d in demands], 
                    capacity[j], 
                    _μ[:,j])
            end
            return (μ, KP_values, 0.0), Inf
        else
            throw(ErrorException("Infeasible subproblem has no dual solution"))
        end
        
    else
        throw(ErrorException("Unexpected dual status"))
    end
end



function calculate_KP_value(costs::Vector{Float64}, demands::Vector{Float64}, capacity::Float64, μ::Vector{Float64})
    n = length(demands)
    
    ratios = Vector{Tuple{Int,Float64}}(undef, n)
    
    negative_count = 0
    for i in 1:n
        ratio = (costs[i]*demands[i] - μ[i])/demands[i]
        if ratio < 0
            negative_count += 1
            ratios[negative_count] = (i, ratio)
        end
    end
    
    resize!(ratios, negative_count)
    sort!(ratios, by=x->x[2])
    
    kp_value = 0.0
    remaining_capacity = capacity
    z = zeros(n)

    for (i, _) in ratios
        if remaining_capacity >= demands[i]
            kp_value += costs[i]*demands[i] - μ[i]
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

function calculate_KP_value_infeasible(costs::Vector{Float64}, demands::Vector{Float64}, capacity::Float64, μ::Vector{Float64})
    n = length(demands)
    
    ratios = Vector{Tuple{Int,Float64}}(undef, n)
    
    negative_count = 0
    for i in 1:n
        ratio =  - μ[i]/demands[i]
        if ratio < 0
            negative_count += 1
            ratios[negative_count] = (i, ratio)
        end
    end
    
    resize!(ratios, negative_count)
    sort!(ratios, by=x->x[2])
    
    kp_value = 0.0
    remaining_capacity = capacity
    z = zeros(n)

    for (i, _) in ratios
        if remaining_capacity >= demands[i]
            kp_value += - μ[i]
            remaining_capacity -= demands[i]
            z[i] = 1.0
        else
            fraction = remaining_capacity / demands[i]
            kp_value +=  - μ[i] * fraction
            z[i] = fraction
            break
        end
    end

    return kp_value
end

# ============================================================================
# Helper functions,for testing, don't remove 
# ============================================================================

# function calculate_KP_value(costs::Vector{Float64}, demands::Vector{Float64}, capacity::Float64, μ::Vector{Float64})
#     # Calculate ratios and store indices
#     n = length(demands)
#     ratios_with_idx = [(i, (costs[i]*demands[i] - μ[i])/demands[i]) for i in 1:n]
#     # Filter out positive ratios and sort remaining by ratio in ascending order
#     negative_ratios = filter(x -> x[2] < 0, ratios_with_idx)
#     sort!(negative_ratios, by=x->x[2], rev=false)
#     # Initialize solution value
#     kp_value = 0.0
#     remaining_capacity = capacity
#     z = zeros(n)
#     # Fill knapsack greedily by best negative ratios
#     for (i, ratio) in negative_ratios
#         if remaining_capacity >= demands[i]
#             # Take whole item
#             kp_value += costs[i]*demands[i] - μ[i]
#             remaining_capacity -= demands[i]
#             z[i] = 1.0
#         else
#             # Take fractional amount
#             fraction = remaining_capacity / demands[i]
#             kp_value += (costs[i]*demands[i] - μ[i]) * fraction
#             z[i] = fraction
#             break
#         end
#     end
#     return kp_value
# end

# function test_knapsack_cut(costs::Vector{Float64}, demands::Vector{Float64}, capacity::Float64, μ::Vector{Float64})
#     model = Model(CPLEX.Optimizer)
#     set_optimizer_attribute(model, MOI.Silent(), true)
    
#     @variable(model, 0 <= z[1:length(demands)] <= 1)
#     @objective(model, Min, sum((demands[j] * costs[j] - μ[j]) * z[j] for j in 1:length(demands)))
#     # @objective(model, Min, sum((μ[j] - demands[j] * costs[j]) * z[j] for j in 1:length(demands)))
#     @constraint(model, sum(demands[j] * z[j] for j in 1:length(demands)) <= capacity)
    
#     optimize!(model)
#     # println(model)
#     z_values = value.(z)
#     z_rounded = [v < 0.001 ? 0.0 : (v > 0.999 ? 1.0 : v) for v in z_values]
    
#     return objective_value(model)
# end


# function test_knapsack_cut_infeasible(costs::Vector{Float64}, demands::Vector{Float64}, capacity::Float64, μ::Vector{Float64})
#     model = Model(CPLEX.Optimizer)
#     # set_optimizer_attribute(model, "InfUnbdInfo", 1)
#     set_optimizer_attribute(model, MOI.Silent(), true)
    
#     @variable(model, 0 <= z[1:length(demands)] <= 1)
#     @objective(model, Min, sum( -μ[j] * z[j] for j in 1:length(demands)))
#     @constraint(model, sum(demands[j] * z[j] for j in 1:length(demands)) <= capacity)
#     # println(model)
#     optimize!(model)
#     # @info termination_status(model)
#     z_values = value.(z)
#     z_rounded = [v < 0.001 ? 0.0 : (v > 0.999 ? 1.0 : v) for v in z_values]
    
#     return objective_value(model)
# end