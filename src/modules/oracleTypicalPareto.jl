export ParetoOracle, ParetoOracleParam, model_reformulation!

mutable struct ParetoOracleParam <: AbstractOracleParam
    rtol::Float64
    atol::Float64
    core_point::Vector{Float64}

    function ParetoOracleParam(; core_point = Float64[], rtol = 1e-9, atol = 1e-9)
        isempty(core_point) && throw(AlgorithmException("Please provide core point"))
        new(rtol, atol, core_point)
    end
end

mutable struct ParetoOracle <: AbstractTypicalOracle
    
    oracle_param::ParetoOracleParam
    solver_param::Dict{String,Any}

    model::Model
    pareto_model::Model

    fixed_x_constraints::Vector{ConstraintRef}

    function ParetoOracle(data::Data;
                            scen_idx::Int=-1, 
                            solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1),
                            oracle_param::ParetoOracleParam = ParetoOracleParam())

        @debug "Building pareto oracle"
        model = Model()
        pareto_model = Model()

        # Define coupling variables and constraints
        @variable(model, x[1:data.dim_x])
        @constraint(model, fix_x, x .== 0)

        assign_attributes!(model, solver_param)
        
        new(oracle_param, solver_param, model, pareto_model, fix_x)
    end

    ParetoOracle() = new()
end

function generate_cuts(oracle::ParetoOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    optimize!(oracle.model)
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    status = dual_status(oracle.model)

    if status == FEASIBLE_POINT
        sub_obj_val = objective_value(oracle.model)

        set_time_limit_sec(oracle.pareto_model, time_limit)
        set_objective_coefficient(oracle.pareto_model, oracle.pareto_model[:z], sub_obj_val)
        set_normalized_coefficient.(oracle.pareto_model[:fix_x], oracle.pareto_model[:z], x_value)
        set_normalized_rhs.(oracle.pareto_model[:fix_x], oracle.oracle_param.core_point)
        optimize!(oracle.pareto_model) 

        a_x = dual.(oracle.pareto_model[:fix_x]) 
        a_t = [-1.0] 
        a_0 = sub_obj_val - a_x'*x_value

        if sub_obj_val >= t_value[1] * (1 + oracle.oracle_param.rtol / tol_normalize) + oracle.oracle_param.atol
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
            if dual_sub_obj_val >= oracle.oracle_param.atol / tol_normalize
                return false, [Hyperplane(a_x, a_t, a_0)], [Inf]
            else
                return true, [Hyperplane(a_x, a_t, a_0)], [Inf]
            end
        end
    else
        throw(UnexpectedModelStatusException("ClassicalOracle: $(status)"))
    end
end

function model_reformulation!(oracle::ParetoOracle) 
    oracle.pareto_model = copy(oracle.model)
    set_optimizer(oracle.pareto_model, Gurobi.Optimizer)

    # Remove original coupling constraint
    if haskey(oracle.pareto_model, :fix_x)
        delete(oracle.pareto_model, oracle.pareto_model[:fix_x])
        unregister(oracle.pareto_model, :fix_x)
    end

    x = oracle.pareto_model[:x]
    @variable(oracle.pareto_model, z)

    # Modify constraints for Unified cut 
    expressions = Dict{Symbol, Vector{Any}}(:inequality => [], :equality => [])
    for (t1, t2) in list_of_constraint_types(oracle.pareto_model)
        t1 == VariableRef && continue
        if t2 == MOI.LessThan{Float64}
            for con in all_constraints(oracle.pareto_model, t1, t2)
                lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                push!(expressions[:inequality], @expression(oracle.pareto_model, - rhs * z + rhs - lhs))
            end
        elseif t2 == MOI.GreaterThan{Float64}
            for con in all_constraints(oracle.pareto_model, t1, t2)
                lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                push!(expressions[:inequality], @expression(oracle.pareto_model, rhs * z + lhs - rhs))
            end
        else
            for con in all_constraints(oracle.pareto_model, t1, t2)
                lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                push!(expressions[:equality], @expression(oracle.pareto_model, rhs * z + lhs - rhs))
            end
        end
        delete.(oracle.pareto_model, all_constraints(oracle.pareto_model, t1, t2))
    end
    
    # Add the constraints, where dual of the equalities will be used for cut generation
    @constraint(oracle.pareto_model, expressions[:inequality] .>= 0)
    @constraint(oracle.pareto_model, a_0, expressions[:equality] .== 0)

    # Change objective function
    @objective(oracle.pareto_model, Min, objective_function(oracle.pareto_model) + z)

    # Constraints to fix master variable
    @constraint(oracle.pareto_model, fix_x, x .+ z .== 0)
    assign_attributes!(oracle.pareto_model, oracle.solver_param)
end