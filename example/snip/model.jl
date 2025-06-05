function update_model!(mip::AbstractMip, data::Data)
    x = mip.model[:x]
    model = mip.model
    
    K = data.problem.num_scenarios
    @variable(model, y[1:data.problem.num_nodes, 1:K] >= 0)
    
    @objective(model, Min, sum(data.problem.scenarios[k][3] * y[data.problem.scenarios[k][1],k] for k in 1:K))

    # Constraints
    # Initial probability at destination nodes
    @constraint(model, [k in 1:K], y[data.problem.scenarios[k][2], k] == 1)
    
    # Probability propagation with/without sensors
    for k in 1:K
        for (idx, (from, to, r, q)) in enumerate(data.problem.D)
            @constraint(model, y[from, k] - q * y[to, k] >= 0)  
            @constraint(model, y[from, k] - r * y[to, k] >= - (r - q) * data.problem.ψ[k][to] * x[idx]) 
        end
        for (from, to, r) in data.problem.A_minus_D
            @constraint(model, y[from, k] - r * y[to, k] >= 0)
        end
    end
    
    # Sensor budget constraint
    @constraint(model, sum(x) <= data.problem.budget)

end


function update_model!(master::AbstractMaster, data::Data)
    x = master.model[:x]

    @constraint(master.model, sum(x) <= data.problem.budget)
end

function update_model!(oracle::AbstractTypicalOracle, data::Data, k::Int)
    model = oracle.model
    x = oracle.model[:x]

    # Variables
    @variable(model, y[1:data.problem.num_nodes] >= 0)
    
    # Objective
    @objective(model, Min, y[data.problem.scenarios[k][1]])
    
    # Initial probability constraints at destination nodes
    @constraint(model, y[data.problem.scenarios[k][2]] == 1)

    # Probability propagation constraints
    # Arcs with potential sensors
    for (idx, (from, to, r, q)) in enumerate(data.problem.D)
        @constraint(model, y[from] - q * y[to] >= 0)
        @constraint(model, y[from] - r * y[to] >= -(r - q) * data.problem.ψ[k][to] * x[idx])
    end
    # Arcs without sensors
    for (from, to, r) in data.problem.A_minus_D
        @constraint(model, y[from] - r * y[to] >= 0)
    end

end

function update_model!(oracle::DisjunctiveOracle, data::Data)
    dcglp = oracle.dcglp 

    @constraint(dcglp, [i=1:2], sum(dcglp[:omega_x][i,:]) <= data.problem.budget * dcglp[:omega_0][i])
end