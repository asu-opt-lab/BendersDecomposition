
function generate_cuts(algo::SequentialBenders, ::StandardCut)
    (coefficients_t, coefficients_x, constant_term), subObjVal = generate_cut_coefficients(algo, StandardCut() )

    cut = @expression(algo.master.model, 
        constant_term + coefficients_x'algo.master.var[:x] + coefficients_t * algo.master.var[:t])

    return cut, subObjVal
end


function generate_cut_coefficients(algo::SequentialBenders, ::StandardCut)
    status = dual_status(algo.sub.model)

    if status == FEASIBLE_POINT
        subObjVal = objective_value(algo.sub.model)

        constant_term = subObjVal
        
        # Coefficients for x variables
        coefficients_x = dual.(algo.sub.fixed_x_constraints)
        
        # Coefficients for t variables
        coefficients_t = -1.0 

        # Adjust constant term
        constant_term -= sum(coefficients_x .* algo.master.x_value)
    
    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(algo.sub.model)
            coefficients_x = dual.(algo.sub.fixed_x_constraints)
            coefficients_t = 0.0
            constant_term = 0.0
        else
            @error "Infeasible subproblem has no dual solution"
            throw(ErrorException("Infeasible subproblem has no dual solution"))
        end
        subObjVal = Inf
    else
        @error "Dual status of subproblem is neither feasible nor infeasible: $status"
        throw(ErrorException("Unexpected dual status"))  
    end

    return (coefficients_t, coefficients_x, constant_term), subObjVal
end
