function generate_cuts(env::BendersEnv, cut_strategy::Union{FatKnapsackCut, SlimKnapsackCut})

    # return vector
    critical_pairs, obj_values, coeff_matrix = generate_cut_coefficients(env.sub, env.master.x_value, cut_strategy)
    
    cuts = Vector{Any}(undef, length(critical_pairs))
    for (index, critical_item) in critical_pairs
        cuts[index] = build_cut(env.master, env.sub, critical_pairs[index], cut_strategy)
    end

    return cuts, obj_values, coeff_matrix
end

# Original function
# function generate_cut_coefficients(sub::AbstractSubProblem, x_value::Vector{Float64}, ::Union{FatKnapsackCut, SlimKnapsackCut})
#     # Input validation
#     J = length(sub.sorted_indices)
   
#     # Pre-allocate arrays for better performance
#     critical_pairs = Vector{Tuple{Int,Int}}(undef, J) # (index of the facility, critical item)
#     obj_values = Vector{Float64}(undef, J)

#     # Process each facility
#     for j in 1:J
#         sorted_indices = sub.sorted_indices[j]
#         c_sorted = sub.sorted_cost_demands[j]
#         x_sorted = x_value[sorted_indices]

#         # Find critical item and calculate contribution
#         k = find_critical_item(c_sorted, x_sorted)
#         critical_pairs[j] = (j, k)

#         # Calculate objective value contribution
#         obj_values[j] = c_sorted[k] - (k > 1 ? sum((c_sorted[k] - c_sorted[i]) * x_sorted[i] for i in 1:k-1) : 0)
#     end

#     return critical_pairs, obj_values
# end

# Testing function for retrieving coefficients
function generate_cut_coefficients(sub::AbstractSubProblem, x_value::Vector{Float64}, ::Union{FatKnapsackCut, SlimKnapsackCut})
    # Input validation
    J = length(sub.sorted_indices)
   
    # Pre-allocate arrays for better performance
    critical_pairs = Vector{Tuple{Int,Int}}(undef, J) # (index of the facility, critical item)
    obj_values = Vector{Float64}(undef, J)

    # Coefficient matrix
    coeff_matrix = zeros(Int64, J ,J)

    # Process each facility
    for j in 1:J
        sorted_indices = sub.sorted_indices[j]
        c_sorted = sub.sorted_cost_demands[j]
        x_sorted = x_value[sorted_indices]

        # Find critical item and calculate contribution
        k = find_critical_item(c_sorted, x_sorted)
        critical_pairs[j] = (j, k)
        sum_obj_value_x = 0
        for i in 1:k-1
            # println("I: ", sorted_indices[i], "J: ", j, "coeff: ",  c_sorted[k] - c_sorted[i])
            sum_obj_value_x += (c_sorted[k] - c_sorted[i]) * x_sorted[i]
            coeff_matrix[sorted_indices[i], j] = -(c_sorted[k] - c_sorted[i])
        end

        # Calculate objective value contribution
        obj_values[j] = c_sorted[k] - sum_obj_value_x
    end
    
    return critical_pairs, obj_values, coeff_matrix
end

function build_cut(master::AbstractMasterProblem, sub::AbstractSubProblem, (index,critical_items)::Tuple{Int,Int}, ::Union{FatKnapsackCut, SlimKnapsackCut})
    k = critical_items
    c_sorted = sub.sorted_cost_demands[index]
    sorted_indices = sub.sorted_indices[index]
    return @expression(master.model, - master.model[:t][index] + c_sorted[k] - sum((c_sorted[k] - c_sorted[i]) * master.model[:x][sorted_indices[i]] for i in 1:k-1))
end

function find_critical_item(c::Vector{Float64}, x::Vector{Float64})
    cumsum_x = cumsum(x)
    k = findfirst(>=(1), cumsum_x)
    return k === nothing ? length(c) : k
end

