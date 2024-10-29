

# Knapsack cuts
function _build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, critical_items::Vector{Int}, ::Union{FatKnapsackCut, SlimKnapsackCut})
    expressions_k = []
    expressions_v = []
    expressions_master = []
    for j in 1:length(critical_items)
        k = critical_items[j]
        c_sorted = sub.sorted_cost_demands[j]
        sorted_indices = sub.sorted_indices[j]
        push!(expressions_k, @expression(dcglp.model, c_sorted[k] * dcglp.model[:k₀] - sum((c_sorted[k] - c_sorted[i]) * dcglp.model[:kₓ][sorted_indices[i]] for i in 1:k-1)))
        push!(expressions_v, @expression(dcglp.model, c_sorted[k] * dcglp.model[:v₀] - sum((c_sorted[k] - c_sorted[i]) * dcglp.model[:vₓ][sorted_indices[i]] for i in 1:k-1)))
        push!(expressions_master, @expression(master.model, c_sorted[k] - sum((c_sorted[k] - c_sorted[i]) * master.model[:x][sorted_indices[i]] for i in 1:k-1)))
    end
    return expressions_k, expressions_v, expressions_master
end

function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, critical_items::Vector{Int}, ::FatKnapsackCut)
    expressions_k, expressions_v, expressions_master = _build_cuts(dcglp, master, sub, critical_items, FatKnapsackCut())
    return [@expression(dcglp.model, expr_k - dcglp.model[:kₜ][j]) for (j, expr_k) in enumerate(expressions_k)], 
           [@expression(dcglp.model, expr_v - dcglp.model[:vₜ][j]) for (j, expr_v) in enumerate(expressions_v)],
           [@expression(master.model, expr_master - master.model[:t][j]) for (j, expr_master) in enumerate(expressions_master)]
end

function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, critical_items::Vector{Int}, ::SlimKnapsackCut)
    expressions_k, expressions_v, expressions_master = _build_cuts(dcglp, master, sub, critical_items, SlimKnapsackCut())
    return [@expression(dcglp.model, expr_k - sum(dcglp.model[:kₜ])) for expr_k in expressions_k],
           [@expression(dcglp.model, expr_v - sum(dcglp.model[:vₜ])) for expr_v in expressions_v],
           [@expression(master.model, expr_master - sum(master.model[:t])) for expr_master in expressions_master]
end