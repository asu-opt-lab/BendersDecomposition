export SNIPMILP, create_milp

struct SNIPMILP <: AbstractMILP
    model::Model
end

function create_milp(data::SNIPData)
    model = Model()
    
    # Set optimizer to silent mode
    set_optimizer_attribute(model, MOI.Silent(), false)

    # Extract dimensions
    K = data.num_scenarios
    A = length(data.D) + length(data.A_minus_D)
    @info "K: $K"
    @info data.num_nodes
    @info "A: $A"
    @info length(data.D)
    # Variables
    @variable(model, y[1:length(data.D)], Bin)  # Binary sensor installation variables
    @variable(model, x[1:data.num_nodes, 1:K] >= 0)  # Probability variables
    
    # Objective: minimize expected probability of undetected traversal
    @objective(model, Min, sum(data.scenarios[k][3] * x[data.scenarios[k][1], k] for k in 1:K))
    
    # Constraints
    # Initial probability at destination nodes
    @constraint(model, [k in 1:K], x[data.scenarios[k][2], k] == 1)
    
    # Probability propagation with/without sensors
    for k in 1:K
        for (idx, (from, to, r, q)) in enumerate(data.D)
            @constraint(model, x[from, k] - q * x[to, k] >= 0)  
            @constraint(model, x[from, k] - r * x[to, k] >= - (r - q) * data.ψ[k][to] * y[idx]) 
        end
        for (from, to, r) in data.A_minus_D
            @constraint(model, x[from, k] - r * x[to, k] >= 0)
        end
    end
    
    # Sensor budget constraint
    @constraint(model, sum(y) <= data.budget)
    
    return SNIPMILP(model)
end
