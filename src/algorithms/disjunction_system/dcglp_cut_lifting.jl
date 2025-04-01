function generate_cuts_lifting(env::BendersEnv, cut_strategy::DisjunctiveCut, zeros_indices, ones_indices)

    sub_obj_val = get_subproblem_value(env.sub, env.master.x_value, cut_strategy.base_cut_strategy) 

    disjunctive_inequality = select_disjunctive_inequality(env.master.x_value)

    update_dcglp_1!(env.dcglp, disjunctive_inequality, cut_strategy)
    
    update_dcglp_lifting!(env.dcglp, zeros_indices, ones_indices)

    solve_dcglp!(env, cut_strategy)
    
    cuts = merge_cuts_lifting(env, cut_strategy, zeros_indices, ones_indices)

    if haskey(env.dcglp.model, :lift_0_k)
        delete.(env.dcglp.model, env.dcglp.model[:lift_0_k])
        unregister(env.dcglp.model, :lift_0_k)
        delete.(env.dcglp.model, env.dcglp.model[:lift_1_k])
        unregister(env.dcglp.model, :lift_1_k)
        delete.(env.dcglp.model, env.dcglp.model[:lift_0_v])
        unregister(env.dcglp.model, :lift_0_v)
        delete.(env.dcglp.model, env.dcglp.model[:lift_1_v])
        unregister(env.dcglp.model, :lift_1_v)
    end

    println("--------------------------------DCGLP cut generation finished--------------------------------")

    return cuts, sub_obj_val
end

function update_dcglp_1!(dcglp::DCGLP, disjunctive_inequality::Tuple{Vector{Int}, Int}, disjunction_system::DisjunctiveCut)
    replace_disjunctive_inequality!(dcglp, disjunctive_inequality, disjunction_system.norm_type)
    update_added_benders_constraints!(dcglp, disjunction_system)
    add_disjunctive_cut!(dcglp)
end

function update_dcglp_lifting!(dcglp::DCGLP, zeros_indices, ones_indices)
    
    # Delete existing lifting constraints if they exist
    if haskey(dcglp.model, :lift_0_k)
        delete.(dcglp.model, dcglp.model[:lift_0_k])
        unregister(dcglp.model, :lift_0_k)
        delete.(dcglp.model, dcglp.model[:lift_1_k])
        unregister(dcglp.model, :lift_1_k)
        delete.(dcglp.model, dcglp.model[:lift_0_v])
        unregister(dcglp.model, :lift_0_v)
        delete.(dcglp.model, dcglp.model[:lift_1_v])
        unregister(dcglp.model, :lift_1_v)
    end
    @constraint(dcglp.model, lift_0_k[i=1:length(zeros_indices)], 0 >= dcglp.model[:kₓ][zeros_indices[i]])
    @constraint(dcglp.model, lift_1_k[i=1:length(ones_indices)], 0 >= dcglp.model[:k₀] - dcglp.model[:kₓ][ones_indices[i]])
    @constraint(dcglp.model, lift_0_v[i=1:length(zeros_indices)], 0 >= dcglp.model[:vₓ][zeros_indices[i]])
    @constraint(dcglp.model, lift_1_v[i=1:length(ones_indices)], 0 >= dcglp.model[:v₀] - dcglp.model[:vₓ][ones_indices[i]])
end

function merge_cuts_lifting(env::BendersEnv, cut_strategy::DisjunctiveCut, zeros_indices, ones_indices)
    γ₀, γₓ, γₜ = generate_disjunctive_cuts(env.dcglp, cut_strategy.cut_strengthening_type, zeros_indices, ones_indices)
    push!(env.dcglp.γ_values, (γ₀, γₓ, γₜ))
    master_disjunctive_cut = @expression(env.master.model, γ₀ + dot(γₓ, env.master.var[:x]) + dot(γₜ, env.master.var[:t]))

    if cut_strategy.include_master_cuts
        push!(env.dcglp.master_cuts, master_disjunctive_cut)
        return env.dcglp.master_cuts
    end
    return master_disjunctive_cut
end

# ============================================================================
# Generate strengthened cuts for the DCGLP
# ============================================================================


function generate_disjunctive_cuts(dcglp::DCGLP, ::PureDisjunctiveCut, zeros_indices, ones_indices)
    # optimize!(dcglp.model)

    γₜ = -dual.(dcglp.γ_constraints[:γₜ])
    γ₀ = dual(dcglp.γ_constraints[:γ₀])
    γₓ = -dual.(dcglp.γ_constraints[:γₓ])

    lift_0_k = dual.(dcglp.model[:lift_0_k])
    lift_1_k = dual.(dcglp.model[:lift_1_k])
    lift_0_v = dual.(dcglp.model[:lift_0_v])
    lift_1_v = dual.(dcglp.model[:lift_1_v])

    # Compute lifted values
    lifted_γ₀ = γ₀ - sum(max.(lift_1_k, lift_1_v))
    lifted_γₓ = copy(γₓ) # Start with γₓ values

    # Update lifted_γₓ for special indices
    for (i, idx) in enumerate(zeros_indices)
        lifted_γₓ[idx] = γₓ[idx] + max(lift_0_k[i],lift_0_v[i])
    end
    for (i, idx) in enumerate(ones_indices)
        lifted_γₓ[idx] = γₓ[idx] - max(lift_1_k[i],lift_1_v[i])
    end
    
    _norm_value = norm(vcat(lifted_γₓ, γₜ), Inf)
    norm_value = max(1.0, _norm_value)
    lifted_γₓ = lifted_γₓ ./ norm_value
    γₜ = γₜ ./ norm_value
    lifted_γ₀ = lifted_γ₀ ./ norm_value
    println("DCGLP Cut Stats: [Norm: $norm_value, Gap: $(1 - 1/norm_value), τ: $(value(dcglp.model[:τ])), τ_lb: $(value(dcglp.model[:τ])/norm_value)]")


    return lifted_γ₀, -lifted_γₓ, -γₜ
end


function generate_disjunctive_cuts(dcglp::DCGLP, ::StrengthenedDisjunctiveCut, zeros_indices, ones_indices)
    # optimize!(dcglp.model)

    γₜ = -dual.(dcglp.γ_constraints[:γₜ])
    γ₀ = dual(dcglp.γ_constraints[:γ₀])
    γₓ = -dual.(dcglp.γ_constraints[:γₓ])

    lift_0_k = dual.(dcglp.model[:lift_0_k])
    lift_1_k = dual.(dcglp.model[:lift_1_k])
    lift_0_v = dual.(dcglp.model[:lift_0_v])
    lift_1_v = dual.(dcglp.model[:lift_1_v])

    # Compute lifted values
    lifted_γ₀ = γ₀ - sum(max.(lift_1_k, lift_1_v))
    lifted_γₓ = copy(γₓ) # Start with γₓ values

    # Update lifted_γₓ for special indices
    for (i, idx) in enumerate(zeros_indices)
        lifted_γₓ[idx] = γₓ[idx] + max(lift_0_k[i],lift_0_v[i])
    end
    for (i, idx) in enumerate(ones_indices)
        lifted_γₓ[idx] = γₓ[idx] - max(lift_1_k[i],lift_1_v[i])
    end

    _norm_value = norm(vcat(lifted_γₓ, γₜ), Inf)
    norm_value = max(1.0, _norm_value)
    lifted_γₓ = lifted_γₓ ./ norm_value
    γₜ = γₜ ./ norm_value
    lifted_γ₀ = lifted_γ₀ ./ norm_value
    println("DCGLP Cut Stats: [Norm: $norm_value, Gap: $(1 - 1/norm_value), τ: $(value(dcglp.model[:τ])), τ_lb: $(value(dcglp.model[:τ])/norm_value)]")


    σ₁::Float64 = dual(dcglp.disjunctive_inequalities_constraints[1]) ./ norm_value
    σ₂::Float64 = dual(dcglp.disjunctive_inequalities_constraints[2]) ./ norm_value
    println("DCGLP Sigma Values: [σ₁: $σ₁, σ₂: $σ₂]")

    if abs(σ₁) <= 1e-6 && abs(σ₂) <= 1e-6
        return lifted_γ₀, -lifted_γₓ, -γₜ
    end

    γ₁ = lifted_γₓ .- dual.(dcglp.model[:conv1]) ./ norm_value
    γ₂ = lifted_γₓ .- dual.(dcglp.model[:conv2]) ./ norm_value
    for (i, idx) in enumerate(zeros_indices)
        γ₁[idx] += (lift_0_k[i] - max(lift_0_k[i],lift_0_v[i])) ./ norm_value
        γ₂[idx] += (lift_0_v[i] - max(lift_0_k[i],lift_0_v[i])) ./ norm_value
    end

    σ_sum = σ₂ + σ₁
    if !iszero(σ_sum)
        m = (γ₁ .- γ₂) / σ_sum
        m_lb = floor.(m)
        m_ub = ceil.(m)
        strengthened_γₓ = min.(γ₁-σ₁*m_lb, γ₂+σ₂*m_ub)
    end
    return lifted_γ₀, -strengthened_γₓ, -γₜ
end