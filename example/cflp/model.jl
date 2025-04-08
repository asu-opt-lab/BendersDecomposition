function update_model!(mip::AbstractMip, data::Data)
    x = mip.model[:x]
    model = mip.model
    
    I, J = data.problem.n_facilities, data.problem.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    
    cost_demands = data.problem.costs .* data.problem.demands'
    @objective(model, Min, data.c_x'* x + sum(cost_demands .* y))
    # Add constraints
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.problem.demands[:] .* y[i,:]) <= data.problem.capacities[i] * x[i])
end


function update_model!(master::AbstractMaster, data::Data)
    x = master.model[:x]

    I = data.problem.n_facilities
    @constraint(master.model, capacity, sum(data.problem.capacities[i] * x[i] for i in 1:I) >= sum(data.problem.demands))
end

function update_model!(oracle::AbstractTypicalOracle, data::Data)
    model = oracle.model
    x = oracle.model[:x]
    other_constraints = oracle.other_constraints

    I, J = data.problem.n_facilities, data.problem.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    # Set objective
    cost_demands = data.problem.costs .* data.problem.demands'
    @objective(model, Min, sum(cost_demands .* y))
    # Add constraints
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.problem.demands[:] .* y[i,:]) <= data.problem.capacities[i] * x[i])

    # Store other constraints
    append!(other_constraints, model[:demand])
    append!(other_constraints, vec(model[:facility_open]))
    append!(other_constraints, model[:capacity])
end
