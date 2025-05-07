function update_model!(mip::AbstractMip, data::Data)
    x = mip.model[:x]
    model = mip.model
    
    I, J = data.problem.n_facilities, data.problem.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    @variable(model, t)
    
    cost_demands = data.problem.costs .* data.problem.demands'
    @objective(model, Min, data.c_x'* x + t)

    # Add constraints
    @constraint(model, t >= sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.problem.demands .* y[i,:]) <= data.problem.capacities[i] * x[i])
    @constraint(model, capacity_total, sum(data.problem.capacities[i] * x[i] for i in 1:I) >= sum(data.problem.demands))
end


function update_model!(master::AbstractMaster, data::Data)
    x = master.model[:x]

    I = data.problem.n_facilities
    @constraint(master.model, capacity, sum(data.problem.capacities[i] * x[i] for i in 1:I) >= sum(data.problem.demands))
end

function update_model!(oracle::AbstractTypicalOracle, data::Data)
    model = oracle.model
    x = oracle.model[:x]

    I, J = data.problem.n_facilities, data.problem.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    # Set objective
    cost_demands = data.problem.costs .* data.problem.demands'
    @objective(model, Min, sum(cost_demands .* y))
    # Add constraints
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.problem.demands[:] .* y[i,:]) <= data.problem.capacities[i] * x[i])
end

function update_model!(oracle::DisjunctiveOracle, data::Data)
    dcglp = oracle.dcglp 

    @constraint(dcglp, [i=1:2], sum(data.problem.capacities[j] * dcglp[:omega_x][i,j] for j in 1:data.problem.n_facilities) >= sum(data.problem.demands) * dcglp[:omega_0][i])
end