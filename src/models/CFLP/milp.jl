export CFLPMILP, create_milp

struct CFLPMILP <: AbstractMILP
    model::Model
end

function create_milp(data::CFLPData)
    model = Model()

    # Set optimizer to silent mode
    set_optimizer_attribute(model, MOI.Silent(), true)

    # Extract problem dimensions
    N, M = data.n_facilities, data.n_customers

    # Define variables
    @variable(model, x[1:N], Bin)
    @variable(model, y[1:N, 1:M] >= 0)

    # Set objective
    @objective(model, Min, 
        sum(data.costs[i,j] * data.demands[j] * y[i,j] for i in 1:N, j in 1:M) + 
        sum(data.fixed_costs[i] * x[i] for i in 1:N)
    )

    # Add constraints
    @constraint(model, demand[j in 1:M], sum(y[:,j]) == 1)
    @constraint(model, facility_open[i in 1:N, j in 1:M], y[i,j] <= x[i])
    @constraint(model, capacity[i in 1:N], sum(data.demands[:] .* y[i,:]) <= data.capacities[i] * x[i])

    return CFLPMILP(model)
end
