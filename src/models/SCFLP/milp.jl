export SCFLPMILP, create_milp

struct SCFLPMILP <: AbstractMILP
    model::Model
end

function create_milp(data::SCFLPData)
    model = Model()

    # Set optimizer to silent mode
    set_optimizer_attribute(model, MOI.Silent(), true)

    # Extract problem dimensions
    N, M, S = data.n_facilities, data.n_customers, data.n_scenarios

    # Define variables
    @variable(model, x[1:N], Bin)
    @variable(model, y[1:N, 1:M, 1:S] >= 0)

    # Set objective
    @objective(model, Min,
        (1/S) * sum(data.costs[i,j] * data.demands[s][j] * y[i,j,s] for i in 1:N, j in 1:M, s in 1:S) +
        # sum(data.costs[i,j] * data.demands[s][j] * y[i,j,s] for i in 1:N, j in 1:M, s in 1:S) +
        sum(data.fixed_costs[i] * x[i] for i in 1:N)
    )

    # Add constraints
    @constraint(model, demand[j in 1:M, s in 1:S], sum(y[:,j,s]) == 1)
    @constraint(model, facility_open[i in 1:N, j in 1:M, s in 1:S], y[i,j,s] <= x[i])
    @constraint(model, capacity[i in 1:N, s in 1:S], sum(data.demands[s][j] * y[i,j,s] for j in 1:M) <= data.capacities[i] * x[i])

    return SCFLPMILP(model)
end
