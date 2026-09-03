function update_master_model!(model::Model, data::SNIPData)
    K = data.num_scenarios
    @variable(model, x[1:length(data.D)], Bin)
    @variable(model, t[1:K] >= -1e6)

    @objective(model, Min, sum(data.scenarios[k][3] * t[k] for k in 1:K))

    @constraint(model, sum(x) <= data.budget)

    return (x = x, ), t
end

function update_sub_model!(model::Model, data::SNIPData; x, scen_idx::Int = 0)
    @variable(model, y[1:data.num_nodes] >= 0)
    
    @objective(model, Min, y[data.scenarios[scen_idx][1]])
    
    @constraint(model, y[data.scenarios[scen_idx][2]] == 1)

    for (idx, (from, to, r, q)) in enumerate(data.D)
        @constraint(model, y[from] - q * y[to] >= 0)
        @constraint(model, y[from] - r * y[to] >= -(r - q) * data.ψ[scen_idx][to] * x[idx])
    end

    for (from, to, r) in data.A_minus_D
        @constraint(model, y[from] - r * y[to] >= 0)
    end
    return nothing
end
