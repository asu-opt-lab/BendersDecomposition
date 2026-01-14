export customize_master_model!, customize_sub_model!, customize_mip_model!

"""
    customize_mip_model!(model::Model, data::SNIPData; optimizer=nothing)

Customize the MIP model for SNIP. If no optimizer is provided, the model will not have an optimizer set.
"""
function customize_mip_model!(model::Model, data::SNIPData; optimizer=nothing)
    if optimizer !== nothing
        set_optimizer(model, optimizer)
    end
    
    K = data.num_scenarios
    @variable(model, x[1:length(data.D)], Bin)
    @variable(model, y[1:data.num_nodes, 1:K] >= 0)
    
    @objective(model, Min, sum(data.scenarios[k][3] * y[data.scenarios[k][1],k] for k in 1:K))

    @constraint(model, [k in 1:K], y[data.scenarios[k][2], k] == 1)
    
    for k in 1:K
        for (idx, (from, to, r, q)) in enumerate(data.D)
            @constraint(model, y[from, k] - q * y[to, k] >= 0)  
            @constraint(model, y[from, k] - r * y[to, k] >= - (r - q) * data.ψ[k][to] * x[idx]) 
        end
        for (from, to, r) in data.A_minus_D
            @constraint(model, y[from, k] - r * y[to, k] >= 0)
        end
    end
    
    @constraint(model, sum(x) <= data.budget)
end

"""
    customize_master_model!(model::Model, data::SNIPData; optimizer=nothing)

Customize the master model for SNIP Benders decomposition.
"""
function customize_master_model!(model::Model, data::SNIPData; optimizer=nothing)
    if optimizer !== nothing
        set_optimizer(model, optimizer)
    end
    
    K = data.num_scenarios
    @variable(model, x[1:length(data.D)], Bin)
    @variable(model, t[1:K] >= -1e6)

    @objective(model, Min, sum(data.scenarios[k][3] * t[k] for k in 1:K))

    @constraint(model, sum(x) <= data.budget)

    return (x = x, ), t
end

"""
    customize_sub_model!(model::Model, data::SNIPData, scen_idx::Int; x, optimizer=nothing)

Customize the subproblem model for SNIP Benders decomposition.
"""
function customize_sub_model!(model::Model, data::SNIPData, scen_idx::Int; x, optimizer=nothing)
    if optimizer !== nothing
        set_optimizer(model, optimizer)
    end

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
end
