export ClassicalOracle

mutable struct ClassicalOracle <: AbstractTypicalOracle
    
    oracle_param::EmptyOracleParam

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}

    function ClassicalOracle(data::Data; 
                             scen_idx::Int=-1, 
                             solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9),
                             oracle_param::EmptyOracleParam = EmptyOracleParam())
        @debug "Building classical oracle"
        model = Model()

        # Define coupling variables and constraints
        @variable(model, x[1:data.dim_x])
        @constraint(model, fix_x, x .== 0)

        assign_attributes!(model, solver_param)
        
        new(oracle_param, model, fix_x)
    end

    ClassicalOracle() = new()
end

function generate_cuts(oracle::ClassicalOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-9, atol = 0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    optimize!(oracle.model)
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end
    
    status = dual_status(oracle.model)
    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)

        a_x = dual.(oracle.fixed_x_constraints) 
        a_t = [-1.0] 
        a_0 = sub_obj_val - a_x'*x_value 
        if sub_obj_val >= t_value[1] * (1 + tol) + atol
            return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        end

    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(oracle.model)
            dual_sub_obj_val = dual_objective_value(oracle.model)
            @info "dual_sub_obj_val = $dual_sub_obj_val"
            a_x = dual.(oracle.fixed_x_constraints)
            a_t = [0.0]
            a_0 = dual_sub_obj_val - a_x'*x_value 
            if dual_sub_obj_val >= atol
                return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
            else
                return true, [Hyperplane(a_x, a_t, a_0)], [Inf]
            end
        end
    else
        throw(UnexpectedModelStatusException("ClassicalOracle: $(status)"))
    end
end




