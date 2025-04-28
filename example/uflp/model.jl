# mip formulation
function update_model!(mip::AbstractMip, data::Data)
    x = mip.model[:x]
    model = mip.model
    
    I, J = data.problem.n_facilities, data.problem.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    @variable(model, t[1:J] >= 0)
    
    cost_demands = data.problem.costs .* data.problem.demands'
    # @objective(model, Min, data.c_x'* x + sum(cost_demands .* y))
    @objective(model, Min, data.c_x'* x + sum(t))
    # Add constraints
    @constraint(model, obj[j in 1:J], t[j] >= sum(cost_demands[:,j] .* y[:,j]))
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
end

function update_model!(master::AbstractMaster, data::Data)
    x = master.model[:x]
    
    @constraint(master.model, sum(x) >= 2)
end

function update_model!(oracle::AbstractTypicalOracle, data::Data)
    model = oracle.model
    x = oracle.model[:x]

    # Extract problem dimensions
    I, J = data.problem.n_facilities, data.problem.n_customers
    
    # Define variables
    @variable(model, y[1:I, 1:J] >= 0)

    # Set objective
    cost_demands = data.problem.costs .* data.problem.demands'
    @objective(model, Min, sum(cost_demands .* y))

    # Add constraints
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open[j in 1:J], y[:, j] .<= x)
end

function update_model!(oracle::DisjunctiveOracle, data::Data)
    dcglp = oracle.dcglp 

    @constraint(dcglp, [i=1:2], sum(dcglp[:omega_x][i,:]) >= 2 * dcglp[:omega_0][i])
end