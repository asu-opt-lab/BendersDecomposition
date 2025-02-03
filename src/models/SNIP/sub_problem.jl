export StandardSNIPSubProblem

abstract type AbstractSNIPSubProblem <: AbstractSubProblem end

mutable struct _StandardSNIPSubProblem <: AbstractSNIPSubProblem
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    other_constraints::Vector{ConstraintRef}
end

"""
    StandardSNIPSubProblem <: AbstractSNIPSubProblem

A mutable struct representing the subproblem for the Sensor Network Installation Problem (SNIP).

# Fields
- `sub_problems::Vector{_StandardSNIPSubProblem}`: A vector of subproblems, one for each scenario
    - Each _StandardSNIPSubProblem contains:
        - `model::Model`: The underlying JuMP optimization model
        - `fixed_x_constraints::Vector{ConstraintRef}`: Constraints fixing the x variables to 0
        - `other_constraints::Vector{ConstraintRef}`: Other constraints for the subproblem

# Related Functions
    create_sub_problem(data::SNIPData, cut_strategy::ClassicalCut)
"""
mutable struct StandardSNIPSubProblem <: AbstractSNIPSubProblem
    sub_problems::Vector{_StandardSNIPSubProblem}
end

function create_sub_problem(data::SNIPData, ::ClassicalCut)

    sub_problems = Vector{_StandardSNIPSubProblem}(undef, data.num_scenarios)
    sub_problems = [_create_sub_problem(data, k) for k in 1:data.num_scenarios]
    
    return StandardSNIPSubProblem(sub_problems)
end

function _create_sub_problem(data::SNIPData, k::Int)
    model = Model()
    
    # Variables
    @variable(model, y[1:data.num_nodes] >= 0)
    @variable(model, x[1:length(data.D)])
    
    # Objective
    @objective(model, Min, y[data.scenarios[k][1]])
    
    # Constraints
    other_constraints = Vector{ConstraintRef}()
    
    # Initial probability constraints at destination nodes
    push!(other_constraints, @constraint(model, y[data.scenarios[k][2]] == 1))

    # Probability propagation constraints
    # Arcs with potential sensors
    for (idx, (from, to, r, q)) in enumerate(data.D)
        push!(other_constraints, @constraint(model, y[from] - q * y[to] >= 0))
        push!(other_constraints, @constraint(model, 
            y[from] - r * y[to] >= -(r - q) * data.ψ[k][to] * x[idx]))
    end
    # Arcs without sensors
    for (from, to, r) in data.A_minus_D
        push!(other_constraints, @constraint(model, y[from] - r * y[to] >= 0))
    end

    # Initial y fixing constraints
    fixed_x_constraints = @constraint(model, x .== 0)
    
    return _StandardSNIPSubProblem(model, fixed_x_constraints, other_constraints)
end