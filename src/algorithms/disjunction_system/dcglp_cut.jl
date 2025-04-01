function generate_cuts(env::BendersEnv, cut_strategy::DisjunctiveCut)

    sub_obj_val = get_subproblem_value(env.sub, env.master.x_value, cut_strategy.base_cut_strategy)
    # sub_obj_val = 1e06 # not use 

    disjunctive_inequality = select_disjunctive_inequality(env.master.x_value)

    update_dcglp!(env.dcglp, disjunctive_inequality, cut_strategy)

    solve_dcglp!(env, cut_strategy)
    
    cuts = merge_cuts(env, cut_strategy)

    println("--------------------------------DCGLP cut generation finished--------------------------------")

    return cuts, sub_obj_val
end

function solve_dcglp!(env::BendersEnv, cut_strategy::DisjunctiveCut)

    log = DCGLPIterationLog()
    state = DCGLPState()

    # x_value, t_value = env.master.x_value, env.master.t_value
    x_value = env.master.x_value
    t_value = objective_value(env.sub.model) #perturbation

    set_optimizer_attribute(env.sub.model, "CPX_PARAM_LPMETHOD", 2)
    set_optimizer_attribute(env.sub.model, "CPX_PARAM_EPOPT", 1e-06)
    set_optimizer_attribute(env.sub.model, "CPX_PARAM_ITLIM", 5000)
    set_optimizer_attribute(env.sub.model, MOI.Silent(), false)

    set_normalized_rhs.(env.dcglp.model[:conx], x_value)
    set_normalized_rhs.(env.dcglp.model[:cont], t_value)
    log.start_time = time()

    while true

        state.iteration += 1

        master_time = @elapsed begin
            k_values, v_values, other_values = solve_and_get_dcglp_values(env.dcglp.model, cut_strategy.norm_type)
        end
        log.master_time += master_time

        sub_k_time = @elapsed begin
            dual_info_k, obj_value_k = generate_cut_coefficients(env.sub, k_values, cut_strategy.base_cut_strategy)
        end
        log.sub_k_time += sub_k_time
        
        sub_v_time = @elapsed begin
            dual_info_v, obj_value_v = generate_cut_coefficients(env.sub, v_values, cut_strategy.base_cut_strategy)
        end
        log.sub_v_time += sub_v_time

        update_bounds!(state, k_values, v_values, other_values, obj_value_k, obj_value_v, t_value, cut_strategy.norm_type)
        
        record_iteration!(log, state)

        cut_strategy.verbose && print_dcglp_iteration_info(state, log)

        is_terminated(state, log) && break
        
        if dual_info_k != []
            add_cuts_k!(env, dual_info_k, cut_strategy)
        end
        if dual_info_v != []
            add_cuts_v!(env, dual_info_v, cut_strategy)
        end
        
    end
    statistics_of_disjunctive_cuts(env)
end


function is_terminated(state, log)
    return state.gap <= 1e-3  || state.UB - state.LB <= 1e-03 || get_total_time(log) >= 200 || state.iteration >= 25 #|| (state.UB_k <= 1e-3 && state.UB_v <= 1e-3) 
end

function print_dcglp_iteration_info(state, log)
    @printf("   Iter: %4d | LB: %8.4f | UB: %8.4f | Gap: %6.2f%% | UB_k: %8.2f | UB_v: %8.2f | Master time: %6.2f | Sub_k time: %6.2f | Sub_v time: %6.2f \n",
           state.iteration, state.LB, state.UB, state.gap, sum(state.UB_k), sum(state.UB_v), log.master_time, log.sub_k_time, log.sub_v_time)
end

# function print_dcglp_iteration_info(state, log)
#     @printf("   Iter: %4d | LB: %12.4f | UB: %11.4f | Gap: %8.2f%% | UB_k: %11.4f | UB_v: %11.4f \n",
#            state.iteration, state.LB, state.UB, state.gap, sum(state.UB_k), sum(state.UB_v))
# end

function update_bounds!(state, k_values, v_values, other_values, obj_value_k, obj_value_v, t_value, norm_type::LNorm)
    state.LB = other_values.τ
    diff_st = obj_value_k .+ obj_value_v .- t_value
    state.UB = update_UB!(state.UB, other_values.sx, diff_st, norm_type)
    state.UB_k = obj_value_k .- k_values.t
    state.UB_v = obj_value_v .- v_values.t
    state.gap = (state.UB - state.LB)/abs(state.UB) * 100
end

update_UB!(UB,_sx,diff_st,::L1Norm) = min(UB,norm([ _sx; diff_st], Inf))
update_UB!(UB,_sx,diff_st,::L2Norm) = min(UB,norm([ _sx; diff_st], 2))
update_UB!(UB,_sx,diff_st,::LInfNorm) = min(UB,norm([ _sx; diff_st], 1))

function merge_cuts(env::BendersEnv, cut_strategy::DisjunctiveCut)
    γ₀, γₓ, γₜ = generate_disjunctive_cuts(env.dcglp, cut_strategy.cut_strengthening_type)
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

function generate_disjunctive_cuts(dcglp::DCGLP, ::PureDisjunctiveCut)
    # optimize!(dcglp.model)
    γₜ = dual.(dcglp.model[:cont])
    γ₀ = dual(dcglp.model[:con0])
    γₓ = dual.(dcglp.model[:conx])
    return γ₀, γₓ, γₜ
end


function generate_disjunctive_cuts(dcglp::DCGLP, ::StrengthenedDisjunctiveCut)
    # optimize!(dcglp.model)
    
    σ₁::Float64 = dual(dcglp.disjunctive_inequalities_constraints[1])
    σ₂::Float64 = dual(dcglp.disjunctive_inequalities_constraints[2])
    γₜ::Float64 = dual(dcglp.γ_constraints[:γₜ])
    γ₀::Float64 = dual(dcglp.γ_constraints[:γ₀])
    γₓ::Vector{Float64} = -dual.(dcglp.γ_constraints[:γₓ])

    println("DCGLP Sigma Values: [σ₁: $σ₁, σ₂: $σ₂]")

    if abs(σ₁) <= 1e-6 && abs(σ₂) <= 1e-6
        return γ₀, -γₓ, γₜ
    end

    γ₁ = γₓ .- dual.(dcglp.model[:conv1])
    γ₂ = γₓ .- dual.(dcglp.model[:conv2])
    σ_sum = σ₂ + σ₁
    if !iszero(σ_sum)
        m = (γ₁ .- γ₂) / σ_sum
        m_lb = floor.(m)
        m_ub = ceil.(m)
        γₓ = min.(γ₁-σ₁*m_lb, γ₂+σ₂*m_ub)
    end
    
    return γ₀, -γₓ, γₜ
end


function statistics_of_disjunctive_cuts(env::BendersEnv)
    coef_dict = Dict{String, Union{Float64, Vector{Float64}}}()
    # system 1
    coef_dict["k₀"] = 0.0
    coef_dict["kₜ"] = 0.0
    coef_dict["kₓ"] = zeros(length(env.dcglp.model[:kₓ]))
    # system 2
    coef_dict["v₀"] = 0.0
    coef_dict["vₜ"] = 0.0
    coef_dict["vₓ"] = zeros(length(env.dcglp.model[:vₓ]))

    println("DCGLP Statistics:")

    # ####### -x>=-1 #######
    # dual_eta1 = dual.(env.dcglp.model[:coneta1])
    # dual_eta2 = dual.(env.dcglp.model[:coneta2])
    # non_zero_indices_dual_eta1 = findall(x -> abs(x) >= 1e-6, dual_eta1)
    # non_zero_indices_dual_eta2 = findall(x -> abs(x) >= 1e-6, dual_eta2)
    # println("Non-zero dual values for eta1 constraints: ", non_zero_indices_dual_eta1, "   dual_eta1: ", dual_eta1[non_zero_indices_dual_eta1])
    # println("Non-zero dual values for eta2 constraints: ", non_zero_indices_dual_eta2, "   dual_eta2: ", dual_eta2[non_zero_indices_dual_eta2])
    # # if non_zero_indices_dual_eta1 != []
    # #     for i in non_zero_indices_dual_eta1
    # #         extract_coef(env.dcglp.model[:coneta1][i], dual_eta1[i], coef_dict)
    # #     end
    # # end
    # # if non_zero_indices_dual_eta2 != []
    # #     for i in non_zero_indices_dual_eta2
    # #         extract_coef(env.dcglp.model[:coneta2][i], dual_eta2[i], coef_dict)
    # #     end
    # # end

    # ####### x>=0 #######
    # dual_conv1 = dual.(env.dcglp.model[:conv1])
    # dual_conv2 = dual.(env.dcglp.model[:conv2])
    # non_zero_indices_dual_conv1 = findall(x -> abs(x) >= 1e-6, dual_conv1)
    # non_zero_indices_dual_conv2 = findall(x -> abs(x) >= 1e-6, dual_conv2)
    # println("Non-zero dual values for conv1 constraints: ", non_zero_indices_dual_conv1)
    # println("Non-zero dual values for conv2 constraints: ", non_zero_indices_dual_conv2)
    # # if non_zero_indices_dual_conv1 != []
    # #     for i in non_zero_indices_dual_conv1
    # #         extract_coef(env.dcglp.model[:conv1][i], dual_conv1[i], coef_dict)
    # #     end
    # # end
    # # if non_zero_indices_dual_conv2 != []
    # #     for i in non_zero_indices_dual_conv2
    # #         extract_coef(env.dcglp.model[:conv2][i], dual_conv2[i], coef_dict)
    # #     end
    # # end

    # ####### k0>=0, v0>=0 #######
    # dual_conk0 = dual.(env.dcglp.model[:conk0])
    # dual_conv0 = dual.(env.dcglp.model[:conv0])
    # println("Dual values for conk0 constraints: ", dual_conk0, "   Dual values for conv0 constraints: ", dual_conv0)
    # # extract_coef(env.dcglp.model[:conk0], dual_conk0, coef_dict)
    # # extract_coef(env.dcglp.model[:conv0], dual_conv0, coef_dict)

    # ####### add_problem_specific_constraints! #######
    # dual_conw1 = dual.(env.dcglp.model[:conw1])
    # dual_conw2 = dual.(env.dcglp.model[:conw2])
    # println("Dual values for conw1 constraints: ", dual_conw1, "   Dual values for conw2 constraints: ", dual_conw2)
    # # extract_coef(env.dcglp.model[:conw1], dual_conw1, coef_dict)
    # # extract_coef(env.dcglp.model[:conw2], dual_conw2, coef_dict)
    
    # ####### disjunctive inequalities #######
    # σ₁::Float64 = dual.(env.dcglp.disjunctive_inequalities_constraints[1])
    # σ₂::Float64 = dual.(env.dcglp.disjunctive_inequalities_constraints[2])
    # println("Dual values for disjunctive inequalities constraints: [σ₁: $σ₁, σ₂: $σ₂]")
    # # extract_coef(env.dcglp.disjunctive_inequalities_constraints[1], σ₁, coef_dict)
    # # extract_coef(env.dcglp.disjunctive_inequalities_constraints[2], σ₂, coef_dict)

    # ####### benders cuts #######
    # for i in 1:length(env.dcglp.dcglp_constraints)
    #     if env.dcglp.dcglp_constraints[i] isa Vector
    #         for j in 1:length(env.dcglp.dcglp_constraints[i])
    #             dual_benders_cuts = dual.(env.dcglp.dcglp_constraints[i][j])
    #             if abs(dual_benders_cuts) >= 1e-6
    #                 println("Dual value : ", dual_benders_cuts)#, "   Cuts: ", env.dcglp.dcglp_constraints[i][j])
    #                 # extract_coef(env.dcglp.dcglp_constraints[i][j], dual_benders_cuts, coef_dict)
    #             end
    #         end
    #     else
    #         dual_benders_cuts = dual(env.dcglp.dcglp_constraints[i])
    #         if abs(dual_benders_cuts) >= 1e-6
    #             println("Dual value : ", dual_benders_cuts)#, "   Cuts: ", env.dcglp.dcglp_constraints[i])
    #             # extract_coef(env.dcglp.dcglp_constraints[i], dual_benders_cuts, coef_dict)
    #         end
    #     end
    # end

    # if haskey(env.dcglp.model, :lift_0_k)
    #     lift_0_k = dual.(env.dcglp.model[:lift_0_k])
    #     lift_1_k = dual.(env.dcglp.model[:lift_1_k])
    #     lift_0_v = dual.(env.dcglp.model[:lift_0_v])
    #     lift_1_v = dual.(env.dcglp.model[:lift_1_v])
    #     non_zero_indices_lift_0_k = findall(x -> abs(x) >= 1e-6, lift_0_k)
    #     non_zero_indices_lift_1_k = findall(x -> abs(x) >= 1e-6, lift_1_k)
    #     non_zero_indices_lift_0_v = findall(x -> abs(x) >= 1e-6, lift_0_v)
    #     non_zero_indices_lift_1_v = findall(x -> abs(x) >= 1e-6, lift_1_v)
    #     println("Non-zero dual values for lift_0_k constraints: ", non_zero_indices_lift_0_k)
    #     println("Non-zero dual values for lift_1_k constraints: ", non_zero_indices_lift_1_k)
    #     println("Non-zero dual values for lift_0_v constraints: ", non_zero_indices_lift_0_v)
    #     println("Non-zero dual values for lift_1_v constraints: ", non_zero_indices_lift_1_v)
    #     # for i in non_zero_indices_lift_0_k
    #     #     extract_coef(env.dcglp.model[:lift_0_k][i], lift_0_k[i], coef_dict)
    #     # end
    #     # for i in non_zero_indices_lift_1_k
    #     #     extract_coef(env.dcglp.model[:lift_1_k][i], lift_1_k[i], coef_dict)
    #     # end
    #     # for i in non_zero_indices_lift_0_v
    #     #     extract_coef(env.dcglp.model[:lift_0_v][i], lift_0_v[i], coef_dict)
    #     # end
    #     # for i in non_zero_indices_lift_1_v
    #     #     extract_coef(env.dcglp.model[:lift_1_v][i], lift_1_v[i], coef_dict)
    #     # end
    # end

    ####### γₜ, γ₀, γₓ #######
    γₜ = dual.(env.dcglp.γ_constraints[:γₜ])
    γ₀ = dual.(env.dcglp.γ_constraints[:γ₀])
    γₓ = dual.(env.dcglp.γ_constraints[:γₓ])
    non_zero_indices_dual_γₓ = findall(x -> abs(x) >= 1e-6, γₓ)
    # println("γₓ: ", γₓ)
    println("Non-zero dual values for γₓ constraints: ", non_zero_indices_dual_γₓ)
    non_zero_indices_dual_γₜ = findall(x -> abs(x) >= 1e-6, γₜ)
    # println("γₜ: ", γₜ)
    println("Non-zero dual values for γₜ constraints: ", non_zero_indices_dual_γₜ)
    disjunctive_cut = @expression(env.master.model, γ₀ + dot(γₓ, env.master.var[:x]) + dot(γₜ, env.master.var[:t]))
    println("disjunctive cut: ", disjunctive_cut)

    # # system 1
    # println("System 1:")
    # non_zero_kₓ = findall(x -> abs(x) >= 1e-6, coef_dict["kₓ"])
    # if isempty(non_zero_kₓ)
    #     println("k₀: ", coef_dict["k₀"], "  kₜ: ", coef_dict["kₜ"], "  kₓ: ", coef_dict["kₓ"])
    # else
    #     println("k₀: ", coef_dict["k₀"], "  kₜ: ", coef_dict["kₜ"], "  Non-zero kₓ: ", [(i, coef_dict["kₓ"][i]) for i in non_zero_kₓ])
    # end
    
    # # system 2
    # println("System 2:")
    # non_zero_vₓ = findall(x -> abs(x) >= 1e-6, coef_dict["vₓ"])
    # if isempty(non_zero_vₓ)
    #     println("v₀: ", coef_dict["v₀"], "  vₜ: ", coef_dict["vₜ"], "  vₓ: ", coef_dict["vₓ"])
    # else
    #     println("v₀: ", coef_dict["v₀"], "  vₜ: ", coef_dict["vₜ"], "  Non-zero vₓ: ", [(i, coef_dict["vₓ"][i]) for i in non_zero_vₓ])
    # end
    
end

function extract_index(var_name)

    m = match(r"\[(\d+)\]", var_name)
    if m !== nothing

        return parse(Int, m.captures[1])
    else
        return nothing
    end
end

function extract_coef(constraint, dual_value, coef_dict)
    constr_object = constraint_object(constraint)
    expr = constr_object.func
    
    if expr isa JuMP.GenericAffExpr
        for (var, coef) in expr.terms
            var_name = name(var)
            base_name = split(var_name, '[')[1]
            index = extract_index(var_name)
            

            if occursin(r"^k0|^k₀", base_name)
                coef_dict["k₀"] += coef * dual_value
                    
            elseif occursin(r"^kt|^kₜ", base_name)
                coef_dict["kₜ"] += coef * dual_value

            elseif occursin(r"^kx|^kₓ", base_name) && index !== nothing
                coef_dict["kₓ"][index] += coef * dual_value
                
            elseif occursin(r"^v0|^v₀", base_name)
                coef_dict["v₀"] += coef * dual_value
                
            elseif occursin(r"^vt|^vₜ", base_name)
                coef_dict["vₜ"] += coef * dual_value
                
            elseif occursin(r"^vx|^vₓ", base_name) && index !== nothing
                coef_dict["vₓ"][index] += coef * dual_value

            else
                error("No match for variable: ", var_name)
            end
        end
    end
end
