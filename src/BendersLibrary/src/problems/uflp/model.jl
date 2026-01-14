export customize_master_model!, customize_sub_model!, customize_mip_model!

"""
    customize_mip_model!(model::Model, data::UFLPData; optimizer=nothing)

Customize the MIP model for UFLP. If no optimizer is provided, the model will not have an optimizer set.
Users must provide their own optimizer using `using CPLEX` or `using Gurobi` and then calling
`get_cplex_optimizer()` or `get_gurobi_optimizer()`.
"""
function customize_mip_model!(model::Model, data::UFLPData; optimizer=nothing)
    
    if optimizer !== nothing
        set_optimizer(model, optimizer)
    end
        
    I, J = data.n_facilities, data.n_customers
    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J] >= 0)
    @variable(model, t[1:J] >= 0)
    
    cost_demands = data.costs .* data.demands'
    @objective(model, Min, data.fixed_costs'* x + sum(t))
    
    @constraint(model, obj[j in 1:J], t[j] >= sum(cost_demands[:,j] .* y[:,j]))
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
end

"""
    customize_master_model!(model::Model, data::UFLPData; optimizer=nothing)

Customize the master model for UFLP Benders decomposition.
"""
function customize_master_model!(model::Model, data::UFLPData; optimizer=nothing)

    if optimizer !== nothing
        set_optimizer(model, optimizer)
    end

    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t >= -1e6)
    
    @constraint(model, sum(x) >= 2)
    
    @objective(model, Min, data.fixed_costs'* x + 1.0 * t)
    
    return (x = x, ), t
end

"""
    customize_sub_model!(model::Model, data::UFLPData, scen_idx::Int; x, optimizer=nothing)

Customize the subproblem model for UFLP Benders decomposition.
"""
function customize_sub_model!(model::Model, data::UFLPData, scen_idx::Int; x, optimizer=nothing) 

    if optimizer !== nothing
        set_optimizer(model, optimizer)
    end

    I, J = data.n_facilities, data.n_customers
    
    @variable(model, y[1:I, 1:J] >= 0)

    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))

    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open[j in 1:J], y[:, j] .<= x)
end
