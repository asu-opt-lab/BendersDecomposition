using JuMP
import BendersX: customize_mip_model!, customize_master_model!, customize_sub_model!

# Requires SCFLPBodurData to be included first.

function customize_mip_model!(model::Model, data::SCFLPBodurData)
    I, J, N = data.n_facilities, data.n_customers, data.n_scenarios

    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J, 1:N] >= 0)

    @objective(model, Min,
        data.fixed_costs' * x +
        (1 / N) * sum(data.costs[i, j] * y[i, j, k] for i in 1:I, j in 1:J, k in 1:N)
    )

    @constraint(model, demand[j in 1:J, k in 1:N],
        sum(y[i, j, k] for i in 1:I) >= data.demands[k][j])

    @constraint(model, capacity[i in 1:I, k in 1:N],
        sum(y[i, j, k] for j in 1:J) <= data.capacities[i] * x[i])

    max_total = maximum(sum(data.demands[k]) for k in 1:N)
    @constraint(model, total_capacity,
        sum(data.capacities[i] * x[i] for i in 1:I) >= max_total)
end

function customize_master_model!(model::Model, data::SCFLPBodurData)
    I, N = data.n_facilities, data.n_scenarios

    @variable(model, x[1:I], Bin)
    @variable(model, t[1:N] >= -1e6)

    @objective(model, Min, data.fixed_costs' * x + fill(1 / N, N)' * t)

    max_total = maximum(sum(data.demands[k]) for k in 1:N)
    @constraint(model, capacity,
        sum(data.capacities[i] * x[i] for i in 1:I) >= max_total)

    return (x = x,), t
end

function customize_sub_model!(model::Model, data::SCFLPBodurData, scen_idx::Int; x)
    I, J = data.n_facilities, data.n_customers

    @variable(model, y[1:I, 1:J] >= 0)

    @objective(model, Min,
        sum(data.costs[i, j] * y[i, j] for i in 1:I, j in 1:J))

    @constraint(model, demand[j in 1:J],
        sum(y[i, j] for i in 1:I) >= data.demands[scen_idx][j])

    @constraint(model, capacity[i in 1:I],
        sum(y[i, j] for j in 1:J) <= data.capacities[i] * x[i])

    return nothing
end
