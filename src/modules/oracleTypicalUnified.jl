export UnifiedOracle, UnifiedOracleParam, model_reformulation!

mutable struct UnifiedOracleParam <: AbstractOracleParam
    rtol::Float64
    atol::Float64

    function UnifiedOracleParam(; rtol = 1e-9, atol = 1e-9)
        new(rtol, atol)
    end
end

mutable struct UnifiedOracle <: AbstractTypicalOracle
    
    oracle_param::UnifiedOracleParam
    solver_param::Dict{String,Any}

    model::Model
    unified_model::Model

    function UnifiedOracle(data::Data;
                            scen_idx::Int=-1, 
                            solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9),
                            oracle_param::UnifiedOracleParam = UnifiedOracleParam())

        @debug "Building unified oracle"
        model = Model()
        unified_model = Model()

        # Define coupling variables and constraints
        @variable(model, x[1:data.dim_x])

        assign_attributes!(model, solver_param)
        
        new(oracle_param, solver_param, model, unified_model)
    end

    UnifiedOracle() = new()
end

function generate_cuts(oracle::UnifiedOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.unified_model, time_limit)

    set_normalized_rhs.(oracle.unified_model[:fix_x_lb], x_value)
    set_normalized_rhs.(oracle.unified_model[:fix_x_ub], -x_value)
    set_normalized_rhs.(oracle.unified_model[:epigraph], -t_value)

    optimize!(oracle.unified_model)

    if termination_status(oracle.unified_model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    end

    @assert termination_status(oracle.unified_model) == OPTIMAL

    a_x = dual.(oracle.unified_model[:fix_x_lb]) .- dual.(oracle.unified_model[:fix_x_ub])
    a_t = [-dual(oracle.unified_model[:epigraph])]
    a_0 = objective_value(oracle.unified_model) - a_x'*x_value + dual(oracle.unified_model[:epigraph])*t_value[1]
    
    return isapprox(dual_objective_value(oracle.unified_model), 0, atol=1e-6) ? (true, [], [NaN]) : (false, [Hyperplane(a_x, a_t, a_0)], [NaN])
end

function model_reformulation!(oracle::UnifiedOracle) 
    oracle.unified_model = copy(oracle.model)

    x = oracle.unified_model[:x]
    @variable(oracle.unified_model, z)

    # Modify constraints for Unified cut 
    expressions = Dict{Symbol, Vector{Any}}(:leq => [], :geq => [], :eq2leq=> [], :eq2geq=> [])
    for (t1, t2) in list_of_constraint_types(oracle.unified_model)
        t1 == VariableRef && continue
        if t2 == MOI.LessThan{Float64}
            for con in all_constraints(oracle.unified_model, t1, t2)
                lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                push!(expressions[:leq], @expression(oracle.unified_model, z .+ rhs .- lhs))
            end
        elseif t2 == MOI.GreaterThan{Float64}
            for con in all_constraints(oracle.unified_model, t1, t2)
                lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                push!(expressions[:geq], @expression(oracle.unified_model, z .+ lhs .- rhs))
            end
        else
            for con in all_constraints(oracle.unified_model, t1, t2)
                lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                push!(expressions[:eq2geq], @expression(oracle.unified_model, z .+ lhs .- rhs))
                push!(expressions[:eq2leq], @expression(oracle.unified_model, z .- lhs .+ rhs))
            end
        end
        delete.(oracle.unified_model, all_constraints(oracle.unified_model, t1, t2))
    end

    # Add the constraints that have a slack
    @constraint(oracle.unified_model, geq, expressions[:geq] .>= 0)
    @constraint(oracle.unified_model, leq, expressions[:leq] .>= 0)
    @constraint(oracle.unified_model, eq2geq, expressions[:eq2geq] .>= 0)
    @constraint(oracle.unified_model, eq2leq, expressions[:eq2leq] .>= 0)

    # Add epigraph constraint
    @constraint(oracle.unified_model, epigraph, z .- objective_function(oracle.unified_model) .>= 0)

    # Change objective function
    @objective(oracle.unified_model, Min, z)

    # Constraints to fix master variable
    @constraint(oracle.unified_model, fix_x_lb, z .+ x .>= 0)
    @constraint(oracle.unified_model, fix_x_ub, z .- x .>= 0)
    
    assign_attributes!(oracle.unified_model, oracle.solver_param)
end