function update_master_model!(model::Model, data::CFLPData)
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t >= -1e6)

    @objective(model, Min, data.fixed_costs'* x + t)

    I = data.n_facilities
    @constraint(model, capacity, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))

    return (x = x, ), t
end

function update_sub_model!(model::Model, data::CFLPData, scen_idx::Int; x)
    I, J = data.n_facilities, data.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    # Set objective
    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))
    # Add constraints
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.demands[:] .* y[i,:]) <= data.capacities[i] * x[i])
    return nothing
end
