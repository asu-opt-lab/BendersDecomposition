export InexactOracle

mutable struct InexactOracle <: AbstractTypicalOracle

    typical_oracle::AbstractOracle
    solver_param::Dict{String,Any}

    function InexactOracle(typical_oracle::AbstractOracle; solver_param::Dict{String,Any} = Dict("solver" => "OSQP"))
        
        new(typical_oracle, solver_param)
    end

    InexactOracle() = new()
end

function generate_cuts(oracle::InexactOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    typical_oracle = oracle.typical_oracle
    model, fixed_x_constraints = typical_oracle.model, typical_oracle.fixed_x_constraints
    typical_solver_param, oracle_param = typical_oracle.solver_param, typical_oracle.oracle_param

    set_time_limit_sec(model, time_limit)
    set_normalized_rhs.(fixed_x_constraints, x_value)
    optimize!(model)
    if termination_status(model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end
    
    if objective_value(model) < t_value[1]
        println("obj of OSQP: $(objective_value(model)), t_value: $(t_value[1])")
        assign_attributes!(model, typical_solver_param)
        is_in_L, hyperplanes, f_x =  generate_cuts(typical_oracle, x_value, t_value)
        assign_attributes!(model, oracle.solver_param)
        return is_in_L, hyperplanes, f_x
    end

    status = dual_status(model)
    
    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(model)
        μ = dual.(model[:demand])
        a_t = [-1.0] 
        
        # Get facility knapsack info
        costs = typical_oracle.facility_knapsack_info.costs
        demands = typical_oracle.facility_knapsack_info.demands
        capacity = typical_oracle.facility_knapsack_info.capacity

        # Calculate KP values for each facility
        KP_values = Vector{Float64}(undef, length(capacity))
        for i in 1:length(capacity)
            KP_values[i] = calculate_KP_value(costs[i,:], demands, capacity[i], μ)
        end

        a_x = KP_values # Vector{Float64}
        a_0 = sum(μ) 

        if sub_obj_val >= t_value[1] * (1 + oracle_param.rtol / tol_normalize)
            # return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
            return false, [Hyperplane(a_x, a_t, a_0)], [NaN]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], deepcopy(t_value)
        end
        
    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(model)
            dual_sub_obj_val = dual_objective_value(model)
            @info "dual_sub_obj_val = $dual_sub_obj_val"
            
            a_x = dual.(fixed_x_constraints)
            a_t = [0.0]
            a_0 = dual_sub_obj_val - a_x' * x_value 
            if dual_sub_obj_val >= oracle_param.atol / tol_normalize
                return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
            else
                return true, [Hyperplane(a_x, a_t, a_0)], [Inf]
            end
        end
    else
        @info("termination_status: $(termination_status(model)), primal: $(primal_status(model)), dual: $status")
        throw(UnexpectedModelStatusException("InexactOracle: $(status)"))
    end
end