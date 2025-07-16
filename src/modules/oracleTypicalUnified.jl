export UnifiedOracle, model_reformulation!

mutable struct UnifiedOracle <: AbstractTypicalOracle
    
    oracle_param::EmptyOracleParam
    solver_param::Dict{String,Any}

    model::Model
    unified_model::Model

    pareto::Bool
    core_point::Vector{Float64}

    function UnifiedOracle(data::Data;
                            scen_idx::Int=-1, 
                            pareto::Bool=false,
                            core_point::Vector{Float64}=zeros(data.dim_x),
                            solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9),
                            oracle_param::EmptyOracleParam = EmptyOracleParam())

        @debug "Building unified oracle"
        model = Model()
        unified_model = Model()

        # Define coupling variables and constraints
        @variable(model, x[1:data.dim_x])
        @variable(model, t[1:1])

        assign_attributes!(model, solver_param)
        
        new(oracle_param, solver_param, model, unified_model, pareto, core_point)
    end

    UnifiedOracle() = new()
end

function generate_cuts(oracle::UnifiedOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.unified_model, time_limit)

    set_normalized_rhs.(oracle.unified_model[:fix_x_lb], x_value)
    set_normalized_rhs.(oracle.unified_model[:fix_t_lb], t_value)
    set_normalized_rhs.(oracle.unified_model[:fix_x_ub], -x_value)
    set_normalized_rhs.(oracle.unified_model[:fix_t_ub], -t_value)

    optimize!(oracle.unified_model)
    
    if termination_status(oracle.unified_model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    @assert termination_status(oracle.unified_model) == OPTIMAL

    fix_x_lb = dual.(oracle.unified_model[:fix_x_lb])
    fix_x_ub = dual.(oracle.unified_model[:fix_x_ub])

    fix_t_lb = dual.(oracle.unified_model[:fix_t_lb])
    fix_t_ub = dual.(oracle.unified_model[:fix_t_ub])

    a_x = fix_x_lb - fix_x_ub
    a_t = fix_t_lb - fix_t_ub
    a_0 = sum(dual.(oracle.unified_model[:a_0]))

    return isapprox(dual_objective_value(oracle.unified_model), 0, atol=1e-6) ? (true, [], [NaN]) : (false, [Hyperplane(a_x, a_t, a_0)], [NaN])
end

function model_reformulation!(oracle::UnifiedOracle) 
    oracle.unified_model = copy(oracle.model)

    x = oracle.unified_model[:x]
    t = oracle.unified_model[:t]
    @variable(oracle.unified_model, z)

    # Modify constraints for Unified cut 
    expressions = Dict{Symbol, Vector{Any}}(:inequality => [], :equality => [])
    for (t1, t2) in list_of_constraint_types(oracle.unified_model)
        t1 == VariableRef && continue
        if t2 == MOI.LessThan{Float64}
            for con in all_constraints(oracle.unified_model, t1, t2)
                lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                push!(expressions[:inequality], @expression(oracle.unified_model, z .+ rhs .- lhs))
            end
        elseif t2 == MOI.GreaterThan{Float64}
            for con in all_constraints(oracle.unified_model, t1, t2)
                lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                push!(expressions[:inequality], @expression(oracle.unified_model, z .+ lhs .- rhs))
            end
        else
            for con in all_constraints(oracle.unified_model, t1, t2)
                lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                push!(expressions[:equality], @expression(oracle.unified_model, z .+ lhs .- rhs))
                push!(expressions[:equality], @expression(oracle.unified_model, z .- lhs .+ rhs))
            end
        end
        delete.(oracle.unified_model, all_constraints(oracle.unified_model, t1, t2))
    end
    
    # Add the constraints, where dual of the equalities will be used for cut generation
    @constraint(oracle.unified_model, expressions[:inequality] .>= 0)
    @constraint(oracle.unified_model, a_0, expressions[:equality] .>= 0)

    # Add epigraph constraint
    @constraint(oracle.unified_model, epigraph, z .+ t .- objective_function(oracle.unified_model) .>= 0)

    # Change objective function
    @objective(oracle.unified_model, Min, z)

    # Constraints to fix master variable
    @constraint(oracle.unified_model, fix_x_lb, z .+ x .>= 0)
    @constraint(oracle.unified_model, fix_x_ub, z .- x .>= 0)
    @constraint(oracle.unified_model, fix_t_lb, z .+ t .>= 0)
    @constraint(oracle.unified_model, fix_t_ub, z .- t .>= 0)

    assign_attributes!(oracle.unified_model, oracle.solver_param)
end