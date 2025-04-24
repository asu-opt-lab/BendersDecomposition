export ClassicalOracle

mutable struct ClassicalOracle <: AbstractTypicalOracle
    
    oracle_param::EmptyOracleParam

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    other_constraints::Vector{ConstraintRef}

    function ClassicalOracle(data::Data; 
                             scen_idx::Int=-1, 
                             solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9),
                             oracle_param::EmptyOracleParam = EmptyOracleParam())
        @debug "Building classical oracle"
        model = Model()

        # Define coupling variables and constraints
        @variable(model, x[1:data.dim_x])
        @constraint(model, fix_x, x .== 0)

        other_constr = Vector{ConstraintRef}()

        assign_attributes!(model, solver_param)
        
        new(oracle_param, model, fix_x, other_constr)
    end

    ClassicalOracle() = new()
end

function generate_cuts(oracle::ClassicalOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol = 1e-8, time_limit = 3600)
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
        if sub_obj_val >= t_value[1] + tol
            return false, [Hyperplane(a_x, a_t, a_0)], [sub_obj_val]
        else
            return true, [Hyperplane(a_x, a_t, a_0)], t_value
        end

    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(oracle.model)
            a_x = dual.(oracle.fixed_x_constraints)
            a_t = [0.0]
            a_0 = dual.(oracle.other_constraints)' * normalized_rhs.(oracle.other_constraints)
            return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
        end
    else
        throw(UnexpectedModelStatusException("ClassicalOracle: $(status)"))
    end
end






