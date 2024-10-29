
function generate_cuts(algo::AbstractBendersAlgorithm, cut_strategy::Union{FatKnapsackCut, SlimKnapsackCut})
    critical_items, obj_values = generate_cut_coefficients(algo.sub, algo.master.x_value, cut_strategy)
    return build_cuts(algo.master, algo.sub, critical_items, cut_strategy), obj_values
end

function generate_cut_coefficients(sub::AbstractSubProblem, x_value::Vector{Float64}, ::Union{FatKnapsackCut, SlimKnapsackCut})
    J = length(sub.sorted_indices)

    obj_values = 0
    critical_items = Int[]
    for j in 1:J
        sorted_indices = sub.sorted_indices[j]
        c_sorted = sub.sorted_cost_demands[j]
        x_sorted = x_value[sorted_indices]
        k = find_critical_item(c_sorted, x_sorted)
        push!(critical_items, k)
        obj_values += c_sorted[k] - (k > 1 ? sum((c_sorted[k] - c_sorted[i]) * x_sorted[i] for i in 1:k-1) : 0)
    end

    return critical_items, obj_values
end

function _build_cuts(master::AbstractMasterProblem, sub::AbstractSubProblem, critical_items::Vector{Int}, ::Union{FatKnapsackCut, SlimKnapsackCut})
    expressions = []
    for j in 1:length(critical_items)
        k = critical_items[j]
        c_sorted = sub.sorted_cost_demands[j]
        sorted_indices = sub.sorted_indices[j]
        push!(expressions, @expression(master.model, c_sorted[k] - sum((c_sorted[k] - c_sorted[i]) * master.model[:x][sorted_indices[i]] for i in 1:k-1)))
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


# Old code for reference
# for j in 1:J
#     c = algo.sub.cost_demands[j]
#     sorted_indices = sortperm(c)
#     c_sorted = c[sorted_indices]
#     x_sorted = cVal[sorted_indices]
#     k = find_critical_item(c_sorted, x_sorted)
    
#     ex = @expression(algo.master.model, 
#         c_sorted[k] - sum((c_sorted[k] - c_sorted[i]) * algo.master.model[:x][sorted_indices[i]] for i in 1:k-1)
#     )
    
#     obj_values += c_sorted[k] - (k > 1 ? sum((c_sorted[k] - c_sorted[i]) * x_sorted[i] for i in 1:k-1) : 0)
#     push!(expressions, ex)
# end