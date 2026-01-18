export UnifiedOracle, UnifiedOracleParam, model_reformulation!

mutable struct UnifiedOracleParam <: AbstractOracleParam
    rtol::Float64
    atol::Float64

    function UnifiedOracleParam(; rtol = 1e-6, atol = 1e-6)
        new(rtol, atol)
    end
end

mutable struct UnifiedOracle <: AbstractTypicalOracle
    
    oracle_param::UnifiedOracleParam

    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    w0::Float64

    function UnifiedOracle(data::Data; 
                            scen_idx::Int=1, 
                            solver_param::Dict{String,Any} = Dict("solver" => "CPLEX", "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, "CPX_PARAM_EPOPT" => 1e-9), 
                            oracle_param::UnifiedOracleParam = UnifiedOracleParam(),
                            w0::Float64 = 1.0)
    
            @debug "Building unified oracle"
            model = Model()

            # Define coupling variables and constraints
            @variable(model, x[1:data.dim_x])

            assign_attributes!(model, solver_param)

            @constraint(model, fix_x, x .== 0)

            new(oracle_param, model, fix_x, w0)
    end

    UnifiedOracle() = new()
end

function generate_cuts(oracle::UnifiedOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)
    
    set_normalized_rhs.(oracle.fixed_x_constraints, x_value)
    set_normalized_rhs.(oracle.model[:epigraph], -t_value)

    optimize!(oracle.model)
    
    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    elseif termination_status(oracle.model) != OPTIMAL
        throw(UnexpectedModelStatusException("UnifiedOracle: $(status). This is likely a numerical issue."))
    end

    a_x = dual.(oracle.fixed_x_constraints)
    a_t = [-dual(oracle.model[:epigraph])]
    a_0 = objective_value(oracle.model) - a_x'*x_value + dual(oracle.model[:epigraph])*t_value[1]
    
    return isapprox(dual_objective_value(oracle.model), 0, atol=oracle.oracle_param.atol) ? (true, [Hyperplane(a_x, a_t, a_0)], [NaN]) : (false, [Hyperplane(a_x, a_t, a_0)], [NaN])
end

function model_reformulation!(model::Model, w0::Float64; x)
    # Weight in dual problem
    @variable(model, z)

    # Modify constraints for Unified cut 
    expressions = Dict{Symbol, Vector{Any}}(:leq => [], :geq => [], :eq2leq=> [], :eq2geq=> [])
    for (t1, t2) in list_of_constraint_types(model)
        if t1 == AffExpr
            if t2 == MOI.LessThan{Float64}
                for con in all_constraints(model, t1, t2)
                    lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                    if any(v -> v in keys(lhs.terms), x)
                        push!(expressions[:leq], @expression(model, z .+ rhs .- lhs))
                    else
                        push!(expressions[:leq], @expression(model, rhs .- lhs))
                    end
                end
            elseif t2 == MOI.GreaterThan{Float64}
                for con in all_constraints(model, t1, t2)
                    lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                    if any(v -> v in keys(lhs.terms), x)
                        push!(expressions[:geq], @expression(model, z .+ lhs .- rhs))
                    else
                        push!(expressions[:geq], @expression(model, lhs .- rhs))
                    end
                end
            elseif t2 == MOI.EqualTo{Float64}
                for con in all_constraints(model, t1, t2)
                    (con in model[:fix_x]) && continue
                    lhs, rhs = JuMP.constraint_object(con).func, normalized_rhs(con)
                    if any(v -> v in keys(lhs.terms), x)
                        push!(expressions[:eq2geq], @expression(model, z .+ lhs .- rhs))
                        push!(expressions[:eq2leq], @expression(model, z .- lhs .+ rhs))
                    else
                        push!(expressions[:eq2geq], @expression(model, lhs .- rhs))
                        push!(expressions[:eq2leq], @expression(model, -lhs .+ rhs))
                    end
                end
            else
                throw(UndefError("An affine constraint is not provided as either ≤, ≥, or =."))
            end
            cons = all_constraints(model, t1, t2)
            to_delete = ConstraintRef[]
            for con in cons
                (con in model[:fix_x]) && continue
                push!(to_delete, con)
            end
            delete.(model, to_delete)
        end
    end

    # Add the constraints that have a slack
    @constraint(model, geq, expressions[:geq] .>= 0)
    @constraint(model, leq, expressions[:leq] .>= 0)
    @constraint(model, eq2geq, expressions[:eq2geq] .>= 0)
    @constraint(model, eq2leq, expressions[:eq2leq] .>= 0)

    # Add epigraph constraint
    @constraint(model, epigraph, w0 * z .- objective_function(model) .>= 0)

    # Change objective function
    @objective(model, Min, z)
end




