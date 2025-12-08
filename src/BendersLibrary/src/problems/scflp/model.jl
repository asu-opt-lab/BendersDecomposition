export customize_master_model!, customize_sub_model!, customize_mip_model!

function customize_mip_model!(model::Model, problem::SCFLPData)
    
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)

    set_optimizer(model, optimizer)
    
    # Extract problem dimensions
    I, J, N = problem.n_facilities, problem.n_customers, problem.n_scenarios
    
    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J, 1:N] >= 0)

    # Set objective
    @objective(model, Min,
        (1/N) * sum(problem.costs[i,j] * problem.demands[s][j] * y[i,j,s] for i in 1:I, j in 1:J, s in 1:N) +
        problem.fixed_costs' * x
    )

    # Add constraints
    @constraint(model, demand[j in 1:J, s in 1:N], sum(y[:,j,s]) == 1)
    @constraint(model, facility_open[i in 1:I, j in 1:J, s in 1:N], y[i,j,s] <= x[i])
    @constraint(model, capacity[i in 1:I, s in 1:N], sum(problem.demands[s][j] * y[i,j,s] for j in 1:J) <= problem.capacities[i] * x[i])
end

function customize_master_model!(model::Model, problem::SCFLPData)
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)

    set_optimizer(model, optimizer)

    I, N = problem.n_facilities, problem.n_scenarios
    @variable(model, x[1:I], Bin)
    @variable(model, t[1:N] >= -1e6)

    @objective(model, Min, problem.fixed_costs'* x + fill(1/N, N)' * t)

    max_demand = maximum(sum(demands) for demands in problem.demands)
    @constraint(model, capacity, sum(problem.capacities[i] * x[i] for i in 1:I) >= max_demand)

    return (x = x, ), t
end

function customize_sub_model!(model::Model, problem::SCFLPData, scen_idx::Int; x) 
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, MOI.Silent() => true)

    set_optimizer(model, optimizer)

    I, J = problem.n_facilities, problem.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    # Set objective
    cost_demands = problem.costs .* problem.demands[scen_idx]'
    @objective(model, Min, sum(cost_demands .* y))
    # Add constraints
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(problem.demands[scen_idx][:] .* y[i,:]) <= problem.capacities[i] * x[i])
end

# """
# To-do
# """
# function update_model!(oracle::DisjunctiveOracle, data::Data{<:SCFLPData})
#     dcglp = oracle.dcglp 

#     I = data.problem.n_facilities
#     max_demand = maximum(sum(demands) for demands in data.problem.demands)
#     @constraint(dcglp, [i=1:2], sum(data.problem.capacities[j] * dcglp[:omega_x][i,j] for j in 1:I) >= max_demand * dcglp[:omega_0][i])
# end
