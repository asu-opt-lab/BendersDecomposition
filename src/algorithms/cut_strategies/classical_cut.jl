function generate_cuts(env::BendersEnv, ::ClassicalCut)
    (coefficients_t, coefficients_x, constant_term), sub_obj_val = generate_cut_coefficients(env.sub, env.master.x_value, ClassicalCut())

    cut = @expression(env.master.model, 
        constant_term + dot(coefficients_x, env.master.var[:x]) + coefficients_t * env.master.var[:t])

    return cut, sub_obj_val
end

function generate_cuts(env::BendersEnv, ::ClassicalCut, scenario::Int)
    (coefficients_t, coefficients_x, constant_term), sub_obj_val = generate_cut_coefficients(env.sub.sub_problems[scenario], env.master.x_value, ClassicalCut())

    cut = @expression(env.master.model, 
        constant_term + dot(coefficients_x, env.master.var[:x]) + coefficients_t * env.master.var[:t][scenario])

    return cut, sub_obj_val
end

function generate_cut_coefficients(sub::AbstractSubProblem, x_value::Vector{Float64}, ::ClassicalCut)
    optimize!(sub.model)
    status = dual_status(sub.model)
    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(sub.model)
        coefficients_x = dual.(sub.fixed_x_constraints)
        coefficients_t = -1.0 
        constant_term = sub_obj_val - dot(coefficients_x, x_value)
        # constant_term = dot(dual.(sub.other_constraints), normalized_rhs.(sub.other_constraints))
        return (coefficients_t, coefficients_x, constant_term), sub_obj_val
        
    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(sub.model)
            coefficients_x = dual.(sub.fixed_x_constraints)
            coefficients_t = 0.0
            constant_term = dot(dual.(sub.other_constraints), normalized_rhs.(sub.other_constraints))
            return (coefficients_t, coefficients_x, constant_term), Inf
        else
            throw(ErrorException("Infeasible subproblem has no dual solution"))
        end
        
    else
        throw(ErrorException("Unexpected dual status"))
    end
end
