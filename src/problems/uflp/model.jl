
function update_mip_model!(model::Model, data::UFLPData)
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

function update_master_model!(model::Model, data::UFLPData)
    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t >= -1e6)

    @constraint(model, sum(x) >= 2)

    @objective(model, Min, data.fixed_costs'* x + 1.0 * t)

    return (x = x, ), t
end

function update_sub_model!(model::Model, data::UFLPData, scen_idx::Int; x)
    I, J = data.n_facilities, data.n_customers

    @variable(model, y[1:I, 1:J] >= 0)

    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))

    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open[j in 1:J], y[:, j] .<= x)
    return nothing
end
