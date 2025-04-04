export generate_cuts

function generate_cuts(env::AbstractBendersEnv, cut_strategy::AbstractCutStrategy)

    dual_info, sub_obj_val = generate_cut_coefficients(env.sub, cut_strategy)
    
    cut = build_expression(env.master, dual_info, cut_strategy) 

    return cut, sub_obj_val
end

function generate_cuts(env::AbstractBendersEnv, cut_strategy::AbstractCutStrategy, scenario::Int)
    
    dual_info, sub_obj_val = generate_cut_coefficients(env.sub.sub_problems[scenario], cut_strategy)
    
    cut = build_expression(env.master, scenario, dual_info, cut_strategy)

    return cut, sub_obj_val
end


########################################################
# Classical Cut
########################################################
struct ClassicalCut <: AbstractCutStrategy end
export ClassicalCut

function generate_cut_coefficients(sub::AbstractSubProblem, ::ClassicalCut)

    status = dual_status(sub.model)
    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(sub.model)
        coef_x = dual.(sub.fixed_x_constraints)
        coef_t = -1.0 
        const_term = sub_obj_val - dot(coef_x, sub.fixed_x_values)
        # constant_term = dot(dual.(sub.other_constraints), normalized_rhs.(sub.other_constraints))
        return (coef_t, coef_x, const_term), sub_obj_val
        
    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(sub.model)
            coef_x = dual.(sub.fixed_x_constraints)
            coef_t = 0.0
            const_term = dot(dual.(sub.other_constraints), normalized_rhs.(sub.other_constraints))
            return (coef_t, coef_x, const_term), Inf
        else
            throw(ErrorException("Infeasible subproblem has no dual solution"))
        end
        
    else
        throw(ErrorException("Unexpected dual status"))
    end
end

function build_expression(master::AbstractMasterProblem, dual_info::Tuple{Float64, Vector{Float64}, Float64}, ::ClassicalCut)
    coef_t, coef_x, const_term = dual_info
    cut = @expression(master.model, 
        const_term + dot(coef_x, master.variables[:integer_variables]) + coef_t * master.variables[:continuous_variables])

    return cut
end

function build_expression(master::AbstractMasterProblem, scenario::Int, dual_info::Tuple{Float64, Vector{Float64}, Float64}, ::ClassicalCut)
    coef_t, coef_x, const_term = dual_info
    cut = @expression(master.model, 
        const_term + dot(coef_x, master.variables[:integer_variables]) + coef_t * master.variables[:continuous_variables][scenario])

    return cut
end

