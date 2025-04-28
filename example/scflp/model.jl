function update_model!(mip::AbstractMip, data::Data)
    x = mip.model[:x]
    model = mip.model
    
    # Extract problem dimensions
    I, J, N = data.problem.n_facilities, data.problem.n_customers, data.problem.n_scenarios
    
    @variable(model, y[1:I, 1:J, 1:N] >= 0)

    # Set objective
    @objective(model, Min,
        (1/N) * sum(data.problem.costs[i,j] * data.problem.demands[s][j] * y[i,j,s] for i in 1:I, j in 1:J, s in 1:N) +
        data.problem.fixed_costs' * x
    )

    # Add constraints
    @constraint(model, demand[j in 1:J, s in 1:N], sum(y[:,j,s]) == 1)
    @constraint(model, facility_open[i in 1:I, j in 1:J, s in 1:N], y[i,j,s] <= x[i])
    @constraint(model, capacity[i in 1:I, s in 1:N], sum(data.problem.demands[s][j] * y[i,j,s] for j in 1:J) <= data.problem.capacities[i] * x[i])
end

function update_model!(master::AbstractMaster, data::Data)
    x = master.model[:x]

    I = data.problem.n_facilities
    max_demand = maximum(sum(demands) for demands in data.problem.demands)
    @constraint(master.model, capacity, sum(data.problem.capacities[i] * x[i] for i in 1:I) >= max_demand)
end

function update_model!(oracle::AbstractTypicalOracle, data::Data, scen_idx::Int)
    model = oracle.model
    x = oracle.model[:x]

    I, J = data.problem.n_facilities, data.problem.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    # Set objective
    cost_demands = data.problem.costs .* data.problem.demands[scen_idx]'
    @objective(model, Min, sum(cost_demands .* y))
    # Add constraints
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.problem.demands[scen_idx][:] .* y[i,:]) <= data.problem.capacities[i] * x[i])
end

function update_model!(oracle::DisjunctiveOracle, data::Data)
    dcglp = oracle.dcglp 

    I = data.problem.n_facilities
    max_demand = maximum(sum(demands) for demands in data.problem.demands)
    @constraint(dcglp, [i=1:2], sum(data.problem.capacities[j] * dcglp[:omega_x][i,j] for j in 1:I) >= max_demand * dcglp[:omega_0][i])
end
