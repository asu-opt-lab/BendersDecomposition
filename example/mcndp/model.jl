function update_model!(mip::AbstractMip, data::Data)
    x = mip.model[:x]
    model = mip.model
    
    N, E, K = data.problem.num_nodes, data.problem.num_arcs, data.problem.num_commodities

    @variable(model, y[1:E, 1:K] >= 0)
    
    @objective(model, Min, data.c_x'* x + sum(data.problem.variable_costs[e]*y[e,k] for e in 1:E, k in 1:K))

    # Add constraints
    for k in 1:K
        origin, destination, q = data.problem.demands[k]
        for node in 1:N
            infl_idx, outfl_idx = findall(t -> t[2] == node, data.problem.arcs), findall(t -> t[1] == node, data.problem.arcs)
            if node == origin
                @constraint(model, sum(y[e,k] for e in outfl_idx) - sum(y[e,k] for e in infl_idx) == q)
            elseif node == destination
                @constraint(model, sum(y[e,k] for e in outfl_idx) - sum(y[e,k] for e in infl_idx) == -q)
            else
                @constraint(model, sum(y[e,k] for e in outfl_idx) - sum(y[e,k] for e in infl_idx) == 0)
            end
        end
    end
    @constraint(model, capacity[e in 1:E], sum(y[e,k] for k in 1:K) <= data.problem.capacities[e] * x[e])
end


function update_model!(master::AbstractMaster, data::Data)
end

function update_model!(oracle::AbstractTypicalOracle, data::Data)
    model = oracle.model
    x = oracle.model[:x]

    N, E, K = data.problem.num_nodes, data.problem.num_arcs, data.problem.num_commodities

    @variable(model, y[1:E, 1:K] >= 0)
    
    @objective(model, Min, sum(data.problem.variable_costs[e]*y[e,k] for e in 1:E, k in 1:K))

    # Add constraints
    for k in 1:K
        origin, destination, q = data.problem.demands[k]
        for node in 1:N
            infl_idx, outfl_idx = findall(t -> t[2] == node, data.problem.arcs), findall(t -> t[1] == node, data.problem.arcs)
            if node == origin
                @constraint(model, sum(y[e,k] for e in outfl_idx) - sum(y[e,k] for e in infl_idx) == q)
            elseif node == destination
                @constraint(model, sum(y[e,k] for e in outfl_idx) - sum(y[e,k] for e in infl_idx) == -q)
            else
                @constraint(model, sum(y[e,k] for e in outfl_idx) - sum(y[e,k] for e in infl_idx) == 0)
            end
        end
    end
    @constraint(model, capacity[e in 1:E], sum(y[e,k] for k in 1:K) <= data.problem.capacities[e] * x[e])
end

function update_model!(oracle::DisjunctiveOracle, data::Data)
end