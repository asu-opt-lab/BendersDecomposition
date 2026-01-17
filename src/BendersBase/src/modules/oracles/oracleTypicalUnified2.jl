export UnifiedOracle, model_reformulation!

const UnifiedOracleParam = BasicOracleParam

mutable struct UnifiedOracle <: AbstractTypicalOracle
    
    param::UnifiedOracleParam

    model::Model
    w0::Float64

    function UnifiedOracle(data::AbstractData, master::Master; 
                            customize = customize_sub_model!,
                            scen_idx::Int=0, 
                            param::UnifiedOracleParam = UnifiedOracleParam(), 
                            w0::Float64 = 1.0)
    
            @debug "Building unified oracle"
            model = Model()

            # Copy the master’s coupling variables into the submodel (with identical axes and symbols)
            x_copy = copy_variables!(model, master.x_tuple)

            # Build the submodel using user-defined customization, passing the copied variables
            customize(model, data, scen_idx; x_copy...)

            # Collect all copied master variables and add linking constraint
            x = var_from_tuple(x_copy)

            # Reformulate subproblem
            model_reformulation!(model, w0; x)

            new(param, model, w0)
    end

    UnifiedOracle() = new()
end

function generate_cuts(oracle::UnifiedOracle, x_value::Vector{Float64}, t_value::Vector{Float64}; tol_normalize = 1.0, time_limit = 3600)
    set_time_limit_sec(oracle.model, time_limit)

    set_normalized_rhs.(oracle.model[:fix_x_lb], x_value)
    set_normalized_rhs.(oracle.model[:fix_x_ub], -x_value)
    set_normalized_rhs.(oracle.model[:epigraph], -t_value)

    optimize!(oracle.model)

    if termination_status(oracle.model) == TIME_LIMIT
        throw(TimeLimitException("Time limit reached during cut generation"))
    elseif termination_status(oracle.model) != OPTIMAL
        throw(UnexpectedModelStatusException("UnifiedOracle: $(termination_status(oracle.model)). This is likely a numerical issue."))
    end
    
    a_x = dual.(oracle.model[:fix_x_lb]) .- dual.(oracle.model[:fix_x_ub])
    a_t = [-dual(oracle.model[:epigraph])]
    a_0 = objective_value(oracle.model) - a_x'*x_value + dual(oracle.model[:epigraph])*t_value[1]

    return isapprox(dual_objective_value(oracle.model), 0, atol=oracle.param.zero_tol) ? (true, [Hyperplane(a_x, a_t, a_0)], [NaN]) : (false, [Hyperplane(a_x, a_t, a_0)], [NaN])
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
                throw(ArgumentError("Unsupported constraint sense: $t2. Expected ≤, ≥, or =."))
            end
            delete.(model, all_constraints(model, t1, t2))
        elseif t1 == VariableRef
            continue
        else
            throw(ArgumentError("Constraint type of $t1 is neither AffExpr nor VariableRef."))
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

    # Linking constraints
    @constraint(model, fix_x_lb, z .+ x .>= 0)
    @constraint(model, fix_x_ub, z .- x .>= 0)
end




