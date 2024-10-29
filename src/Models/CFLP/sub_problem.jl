abstract type AbstractCFLPSubProblem end

struct CFLPStandardSubProblem <: AbstractCFLPSubProblem
    model::Model
    constraints::Dict{Symbol, Vector{ConstraintRef}}
    variables::Dict{Symbol, Vector{VariableRef}}
    objective::AffExpr
    data::Any
    solution::Dict{Symbol, Any}
end

struct CFLPNormalizedSubProblem <: AbstractCFLPSubProblem
    model::Model
    constraints::Dict{Symbol, Vector{ConstraintRef}}
    variables::Dict{Symbol, Vector{VariableRef}}
    objective::AffExpr
    data::Any
    solution::Dict{Symbol, Any}
end

struct CFLPNoModelSubProblem <: AbstractCFLPSubProblem
    data::Any
    solution::Dict{Symbol, Any}
end

function create_sub_problem(::Type{T}, data::CFLPData, args...; kwargs...) where T <: AbstractCFLPSubProblem
    if T == CFLPStandardSubProblem
        return create_standard_sub_problem(data, args...; kwargs...)
    elseif T == CFLPNormalizedSubProblem
        return create_normalized_sub_problem(data, args...; kwargs...)
    elseif T == CFLPNoModelSubProblem
        return create_no_model_sub_problem(data, args...; kwargs...)
    else
        error("Unknown sub-problem type: $T")
    end
end

function create_standard_sub_problem(data::CFLPData, x::Vector{Float64}, solver::Symbol=:Gurobi)
    # Create the model with the specified solver
    model = if solver == :CPLEX
        Model(CPLEX.Optimizer)
    elseif solver == :Gurobi
        Model(Gurobi.Optimizer)
    else
        error("Unsupported solver: $solver")
    end
    
    # Set the optimizer to silent mode
    set_optimizer_attribute(model, MOI.Silent(), true)
    
    println("########### Building Standard Subproblem ###########")

    # Extract problem dimensions
    N, M = data.n_facilities, data.n_customers
    
    # Define variables
    @variable(model, y[1:N, 1:M] >= 0)
    @variable(model, x[1:N])

    # Set objective
    @objective(model, Min, sum(data.costs[i,j] * data.demands[j] * y[i,j] for i in 1:N, j in 1:M))

    # Add constraints
    @constraint(model, demand_satisfaction[j in 1:M], sum(y[i,j] for i in 1:N) == 1)
    @constraint(model, capacity[i in 1:N], sum(data.demands[j] * y[i,j] for j in 1:M) <= data.capacities[i] * x[i])
    @constraint(model, facility_open[i in 1:N, j in 1:M], y[i,j] <= x[i])

    # Store all affine constraints
    all_constraints = [
        all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.GreaterThan{Float64});
        all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64});
        all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.EqualTo{Float64});
    ]

    # Initialize constraint and RHS vectors
    constraints = ConstraintRef[]
    rhs = Float64[]

    # Populate constraint and RHS vectors
    for c in all_constraints
        push!(constraints, c)
        push!(rhs, normalized_rhs(c))
    end

    # Add constraints to fix x variables
    fixed_x_constraints = @constraint(model, [i in 1:N], x[i] == 0)

    # Create and return the CFLPStandardSubProblem struct
    return CFLPStandardSubProblem(
        model,
        Dict(:all => constraints, :fixed_x => fixed_x_constraints),
        Dict(:y => y, :x => x),
        objective_function(model),
        data,
        Dict{Symbol, Any}()
    )
end

function create_normalized_sub_problem(data::CFLPData, solver::Symbol=:Gurobi)
    # Implementation for normalized sub-problem
    # ...
end

function create_no_model_sub_problem(data::CFLPData, args...; kwargs...)
    # Implementation for sub-problem without modeling
    return CFLPNoModelSubProblem(data, Dict{Symbol, Any}())
end
