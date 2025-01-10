export StandardSNIPSubProblem

abstract type AbstractSNIPSubProblem <: AbstractSubProblem end

mutable struct StandardSNIPSubProblem <: AbstractSNIPSubProblem
    model::Model
    fixed_x_constraints::Vector{ConstraintRef}
    other_constraints::Vector{ConstraintRef}
end

function create_sub_problem(data::SNIPData, ::ClassicalCut)
    model = Model()
    
    # Variables
    @variable(model, y[1:data.num_nodes, 1:data.num_scenarios] >= 0)
    @variable(model, x[1:length(data.D)])
    
    # Objective
    @objective(model, Min, sum(data.scenarios[k][3] * y[data.scenarios[k][1], k] for k in 1:data.num_scenarios))
    
    # Constraints
    other_constraints = Vector{ConstraintRef}()
    
    # Initial probability constraints at destination nodes
    for k in 1:data.num_scenarios
        push!(other_constraints, @constraint(model, y[data.scenarios[k][2], k] == 1))
    end
    
    # Probability propagation constraints
    for k in 1:data.num_scenarios
        # Arcs with potential sensors
        for (idx, (from, to, r, q)) in enumerate(data.D)
            push!(other_constraints, @constraint(model, y[from, k] - q * y[to, k] >= 0))
            push!(other_constraints, @constraint(model, 
                y[from, k] - r * y[to, k] >= -(r - q) * data.ψ[k][to] * x[idx]))
        end
        # Arcs without sensors
        for (from, to, r) in data.A_minus_D
            push!(other_constraints, @constraint(model, y[from, k] - r * y[to, k] >= 0))
        end
    end
    
    # Initial y fixing constraints
    fixed_x_constraints = @constraint(model, x .== 0)
    
    return StandardSNIPSubProblem(model, fixed_x_constraints, other_constraints)
end 