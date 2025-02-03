function solve_and_get_dcglp_values(model::Model, norm_type::LNorm)
    optimize!(model)
    k_values = (constant = value(model[:k₀]), x = value.(model[:kₓ]), t = value.(model[:kₜ]))
    v_values = (constant = value(model[:v₀]), x = value.(model[:vₓ]), t = value.(model[:vₜ]))
    other_values = (τ = value(model[:τ]), sx = value.(model[:sx]))
    return k_values, v_values, other_values
end

function generate_cut_coefficients(sub::AbstractSubProblem, korv_values, base_cut_strategy::CutStrategy)
    if isapprox(korv_values.constant, 0.0, atol=1e-05)
        return  [], max.(0, korv_values.t)
    end
    input_x = @. abs(korv_values.x / korv_values.constant) # abs value
    solve_sub!(sub, input_x)
    dual_values, obj_value = generate_cut_coefficients(sub, input_x, base_cut_strategy)
    obj_value *= korv_values.constant
    dual_values, obj_value = correct_cut_and_obj_values!(dual_values, obj_value, korv_values.t)
    return  dual_values, obj_value
end

function correct_cut_and_obj_values!(dual_values::Any, obj_value::Float64, t_values::Float64)
    if obj_value <= t_values + 1e-04
        dual_values = []
        obj_value = t_values
    end
    return dual_values, obj_value
end

function correct_cut_and_obj_values!(dual_values::Any, obj_value::Vector{Float64}, t_values::Vector{Float64})
    valid_indices = trues(length(t_values))

    for i in eachindex(t_values)
        if obj_value[i] <= t_values[i] + 1e-04
            valid_indices[i] = false
            obj_value[i] = t_values[i]
        end
    end
    
    if !any(valid_indices)
        dual_values = []
    else
        deleteat!(dual_values, findall(.!valid_indices))
    end
    return dual_values, obj_value
end

# ============================================================================
# Different types of cuts for the DCGLP
# ============================================================================

function add_cuts_k!(env::BendersEnv, dual_info_k, cut_strategy::CutStrategy)
    cuts_k, cuts_v, cuts_master = build_cuts(env.dcglp, env.master, env.sub, dual_info_k, cut_strategy.base_cut_strategy)
    push!(env.dcglp.dcglp_constraints, @constraint(env.dcglp.model, 0 .>= cuts_k))
    if cut_strategy.use_two_sided_cuts 
        push!(env.dcglp.dcglp_constraints, @constraint(env.dcglp.model, 0 .>= cuts_v))
    end
    if cut_strategy.include_master_cuts 
        for cut in cuts_master
            push!(env.dcglp.master_cuts, cut)
        end
    end
end

function add_cuts_v!(env::BendersEnv, dual_info_v, cut_strategy::CutStrategy)
    cuts_k, cuts_v, cuts_master = build_cuts(env.dcglp, env.master, env.sub, dual_info_v, cut_strategy.base_cut_strategy)
    push!(env.dcglp.dcglp_constraints, @constraint(env.dcglp.model, 0 .>= cuts_v))
    if cut_strategy.use_two_sided_cuts 
        push!(env.dcglp.dcglp_constraints, @constraint(env.dcglp.model, 0 .>= cuts_k))
    end
    if cut_strategy.include_master_cuts 
        for cut in cuts_master
            push!(env.dcglp.master_cuts, cut)
        end
    end
end

# knapsack cuts
# optimality cuts
function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, coeff_info::Tuple{Vector{Float64}, Vector{Float64}, Float64}, ::KnapsackCut)
    μ, KP_values, coeff_t = coeff_info
    cut_k = @expression(dcglp.model, coeff_t * dcglp.model[:kₜ] + μ_term(sub, μ) * dcglp.model[:k₀] + sum(KP_values .* dcglp.model[:kₓ]))
    cut_v = @expression(dcglp.model, coeff_t * dcglp.model[:vₜ] + μ_term(sub, μ) * dcglp.model[:v₀] + sum(KP_values .* dcglp.model[:vₓ]))
    cut = @expression(master.model, coeff_t * master.var[:t] + μ_term(sub, μ) + sum(KP_values .* master.var[:x]))
    return cut_k, cut_v, cut
end

# feasibility cuts
function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, coeff_info::Tuple{Float64, Vector{Float64}, Float64}, ::KnapsackCut)
    coefficients_t, coefficients_x, constant_term = coeff_info
    cut_k = @expression(dcglp.model, constant_term*dcglp.model[:k₀] + sum(coefficients_x .* dcglp.model[:kₓ]) + coefficients_t * dcglp.model[:kₜ])
    cut_v = @expression(dcglp.model, constant_term*dcglp.model[:v₀] + sum(coefficients_x .* dcglp.model[:vₓ]) + coefficients_t * dcglp.model[:vₜ])
    cut = @expression(master.model, constant_term + sum(coefficients_x .* master.model[:x]) + coefficients_t * master.model[:t])
    return cut_k, cut_v, cut
end


# fat or slim Knapsack cuts
function _build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, (index, critical_items)::Tuple{Int,Int}, ::Union{FatKnapsackCut, SlimKnapsackCut})
    k = critical_items
    c_sorted = sub.sorted_cost_demands[index]
    sorted_indices = sub.sorted_indices[index]
    expressions_k = @expression(dcglp.model, -dcglp.model[:kₜ][index] + c_sorted[k] * dcglp.model[:k₀] - sum((c_sorted[k] - c_sorted[i]) * dcglp.model[:kₓ][sorted_indices[i]] for i in 1:k-1))
    expressions_v = @expression(dcglp.model, -dcglp.model[:vₜ][index] + c_sorted[k] * dcglp.model[:v₀] - sum((c_sorted[k] - c_sorted[i]) * dcglp.model[:vₓ][sorted_indices[i]] for i in 1:k-1))
    expressions_master = @expression(master.model, -master.model[:t][index] + c_sorted[k] - sum((c_sorted[k] - c_sorted[i]) * master.model[:x][sorted_indices[i]] for i in 1:k-1))
    return expressions_k, expressions_v, expressions_master
end

function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, critical_pairs::Vector{Tuple{Int,Int}}, ::FatKnapsackCut)
    expressions_k = Vector{Any}(undef, length(critical_pairs))
    expressions_v = Vector{Any}(undef, length(critical_pairs))
    expressions_master = Vector{Any}(undef, length(critical_pairs))
    for (i,(index, critical_item)) in enumerate(critical_pairs)
        expressions_k[i], expressions_v[i], expressions_master[i] = _build_cuts(dcglp, master, sub, (index, critical_item), FatKnapsackCut())
    end
    return expressions_k, expressions_v, expressions_master
end

# classical cut
function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, coeff_info::Tuple{Float64, Vector{Float64}, Float64}, ::ClassicalCut)
    coefficients_t, coefficients_x, constant_term = coeff_info
    cut_k = @expression(dcglp.model, constant_term*dcglp.model[:k₀] + sum(coefficients_x .* dcglp.model[:kₓ]) + coefficients_t * dcglp.model[:kₜ])
    cut_v = @expression(dcglp.model, constant_term*dcglp.model[:v₀] + sum(coefficients_x .* dcglp.model[:vₓ]) + coefficients_t * dcglp.model[:vₜ])
    cut = @expression(master.model, constant_term + sum(coefficients_x .* master.model[:x]) + coefficients_t * master.model[:t])
    return cut_k, cut_v, cut
end