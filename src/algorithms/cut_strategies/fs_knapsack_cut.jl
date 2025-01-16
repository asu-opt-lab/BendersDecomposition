function generate_cuts(env::BendersEnv, cut_strategy::Union{FatKnapsackCut, SlimKnapsackCut})
    
    # Generate cut coefficients
    critical_items, obj_values = generate_cut_coefficients(env.sub, env.master.x_value, cut_strategy)
    
    # Build and return cuts
    cuts = build_cuts(env.master, env.sub, critical_items, cut_strategy)
    return cuts, obj_values
end

function generate_cut_coefficients(sub::AbstractSubProblem, x_value::Vector{Float64}, 
                                 ::Union{FatKnapsackCut, SlimKnapsackCut})
    # Input validation
    J = length(sub.sorted_indices)
   
    # Pre-allocate arrays for better performance
    critical_pairs = Vector{Tuple{Int,Int}}(undef, J)
    obj_values = Vector{Float64}(undef, J)

    # Process each facility
    for j in 1:J
        sorted_indices = sub.sorted_indices[j]
        c_sorted = sub.sorted_cost_demands[j]
        x_sorted = x_value[sorted_indices]

        # Find critical item and calculate contribution
        k = find_critical_item(c_sorted, x_sorted)
        critical_pairs[j] = (j, k)

        # Calculate objective value contribution
        obj_values[j] = c_sorted[k] - (k > 1 ? sum((c_sorted[k] - c_sorted[i]) * x_sorted[i] for i in 1:k-1) : 0)
    end

    return critical_items, obj_values
end

function _build_cuts(master::AbstractMasterProblem, sub::AbstractSubProblem, critical_items::Vector{Int}, ::Union{FatKnapsackCut, SlimKnapsackCut})
    expressions = Vector{Any}(undef, length(critical_items))
    for j in 1:length(critical_items)
        k = critical_items[j]
        c_sorted = sub.sorted_cost_demands[j]
        sorted_indices = sub.sorted_indices[j]
        expressions[j] = @expression(master.model, c_sorted[k] - sum((c_sorted[k] - c_sorted[i]) * master.model[:x][sorted_indices[i]] for i in 1:k-1))
    end
    return expressions
end

function build_cuts(master::AbstractMasterProblem, sub::AbstractSubProblem, critical_items::Vector{Int}, ::FatKnapsackCut)
    expressions = _build_cuts(master, sub, critical_items, FatKnapsackCut())
    return [@expression(master.model, expr - master.model[:t][j]) for (j, expr) in enumerate(expressions)]
end

function build_cuts(master::AbstractMasterProblem, sub::AbstractSubProblem, critical_items::Vector{Int}, ::SlimKnapsackCut)
    expressions = _build_cuts(master, sub, critical_items, SlimKnapsackCut())
    return [@expression(master.model, sum(expressions) - sum(master.model[:t]))]
end


function find_critical_item(c::Vector{Float64}, x::Vector{Float64})
    cumsum_x = cumsum(x)
    k = findfirst(>=(1), cumsum_x)
    return k === nothing ? length(c) : k
end

