function solve_and_get_dcglp_values(model::Model, norm_type::LNorm)
    optimize!(model)
    # GRBprintquality(model)
    # GRBgetenv(model)
    # kappa = MOI.get(model, Gurobi.ModelAttribute("KappaExact"))
    # @info "KappaExact: $kappa"
    k_values = (constant = value(model[:k₀]), x = value.(model[:kₓ]), t = value.(model[:kₜ]))
    v_values = (constant = value(model[:v₀]), x = value.(model[:vₓ]), t = value.(model[:vₜ]))
    other_values = (τ = value(model[:τ]), sx = value.(model[:sx]))

    # @info "k_values.x ./ k_values.constant: $(k_values.x ./ k_values.constant)"
    # @info "k_values.t ./ k_values.constant: $(k_values.t ./ k_values.constant)"
    # @info "v_values.x ./ v_values.constant: $(v_values.x ./ v_values.constant)"
    # @info "v_values.t ./ v_values.constant: $(v_values.t ./ v_values.constant)"
    # @info "k_values.t + v_values.t: $(k_values.t + v_values.t)"
    # @info "k_values.x + v_values.x: $(k_values.x + v_values.x)"
    return k_values, v_values, other_values
end

function generate_cut_coefficients(sub::AbstractSubProblem, korv_values, base_cut_strategy::CutStrategy)
    if isapprox(korv_values.constant, 0.0, atol=1e-05)
        return false, [], korv_values.t
    end
    input = @. abs(korv_values.x / korv_values.constant) # abs value
    solve_sub!(sub, input)
    dual_values, obj_value = generate_cut_coefficients(sub, input, base_cut_strategy)
    obj_value *= korv_values.constant
    if obj_value <= korv_values.t - 1e-04
        return false, [], korv_values.t
    end
    return true, dual_values, obj_value
end

function generate_cut_coefficients(sub::AbstractSubProblem, korv_values::NamedTuple, base_cut_strategy::FatKnapsackCut)
    if isapprox(korv_values.constant, 0.0, atol=1e-05)
        return false, [], korv_values.t
    end
    input = @. abs(korv_values.x / korv_values.constant) # abs value
    solve_sub!(sub, input)
    dual_values, obj_value = generate_cut_coefficients(sub, input, base_cut_strategy)
    obj_value *= korv_values.constant
    if sum(obj_value) <= sum(korv_values.t) - 1e-04
        return false, [], korv_values.t
    end
    return true, dual_values, obj_value
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
        push!(env.dcglp.master_cuts, cuts_master)
    end
end

function add_cuts_v!(env::BendersEnv, dual_info_v, cut_strategy::CutStrategy)
    cuts_k, cuts_v, cuts_master = build_cuts(env.dcglp, env.master, env.sub, dual_info_v, cut_strategy.base_cut_strategy)
    push!(env.dcglp.dcglp_constraints, @constraint(env.dcglp.model, 0 .>= cuts_v))
    if cut_strategy.use_two_sided_cuts 
        push!(env.dcglp.dcglp_constraints, @constraint(env.dcglp.model, 0 .>= cuts_k))
    end
    if cut_strategy.include_master_cuts 
        push!(env.dcglp.master_cuts, cuts_master)
    end
end

# knapsack cuts

function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, coeff_info::Tuple{Vector{Float64}, Vector{Float64}, Float64}, ::KnapsackCut)
    μ, KP_values, coeff_t = coeff_info
    cut_k = @expression(dcglp.model, coeff_t * dcglp.model[:kₜ] + sum(μ) * dcglp.model[:k₀] + sum(KP_values .* dcglp.model[:kₓ]))
    cut_v = @expression(dcglp.model, coeff_t * dcglp.model[:vₜ] + sum(μ) * dcglp.model[:v₀] + sum(KP_values .* dcglp.model[:vₓ]))
    cut = @expression(master.model, coeff_t * master.var[:t] + sum(μ) + sum(KP_values .* master.var[:x]))
    return [cut_k], [cut_v], [cut]
end

function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, coeff_info::Tuple{Float64, Vector{Float64}, Float64}, ::KnapsackCut)
    coefficients_t, coefficients_x, constant_term = coeff_info
    cut_k = @expression(dcglp.model, constant_term*dcglp.model[:k₀] + sum(coefficients_x .* dcglp.model[:kₓ]) + coefficients_t * dcglp.model[:kₜ])
    cut_v = @expression(dcglp.model, constant_term*dcglp.model[:v₀] + sum(coefficients_x .* dcglp.model[:vₓ]) + coefficients_t * dcglp.model[:vₜ])
    cut = @expression(master.model, constant_term + sum(coefficients_x .* master.model[:x]) + coefficients_t * master.model[:t])
    return [cut_k], [cut_v], [cut]
end


# fat or slim Knapsack cuts
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
    return [@expression(dcglp.model, sum(expressions_k) - sum(dcglp.model[:kₜ]))],
           [@expression(dcglp.model, sum(expressions_v) - sum(dcglp.model[:vₜ]))],
           [@expression(master.model, sum(expressions_master) - sum(master.model[:t]))]
end

# standard cut
function build_cuts(dcglp::AbstractDCGLP, master::AbstractMasterProblem, sub::AbstractSubProblem, coeff_info::Tuple{Float64, Vector{Float64}, Float64}, ::ClassicalCut)
    coefficients_t, coefficients_x, constant_term = coeff_info
    cut_k = @expression(dcglp.model, constant_term*dcglp.model[:k₀] + sum(coefficients_x .* dcglp.model[:kₓ]) + coefficients_t * dcglp.model[:kₜ])
    cut_v = @expression(dcglp.model, constant_term*dcglp.model[:v₀] + sum(coefficients_x .* dcglp.model[:vₓ]) + coefficients_t * dcglp.model[:vₜ])
    cut = @expression(master.model, constant_term + sum(coefficients_x .* master.model[:x]) + coefficients_t * master.model[:t])
    return [cut_k], [cut_v], [cut]
end