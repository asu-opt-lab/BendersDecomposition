export customize_master_model!, customize_sub_model!, customize_mip_model!

function customize_mip_model!(model::Model, data::SCFLPData)
    
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)

    set_optimizer(model, optimizer)
    
    # Extract dimensions
    I, J, N = data.n_facilities, data.n_customers, data.n_scenarios
    
    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J, 1:N] >= 0)

    # Set objective
    @objective(model, Min,
        (1/N) * sum(data.costs[i,j] * data.demands[s][j] * y[i,j,s] for i in 1:I, j in 1:J, s in 1:N) +
        data.fixed_costs' * x
    )

    # Add constraints
    @constraint(model, demand[j in 1:J, s in 1:N], sum(y[:,j,s]) == 1)
    @constraint(model, facility_open[i in 1:I, j in 1:J, s in 1:N], y[i,j,s] <= x[i])
    @constraint(model, capacity[i in 1:I, s in 1:N], sum(data.demands[s][j] * y[i,j,s] for j in 1:J) <= data.capacities[i] * x[i])
end

function customize_master_model!(model::Model, data::SCFLPData)
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)

    set_optimizer(model, optimizer)

    I, N = data.n_facilities, data.n_scenarios
    @variable(model, x[1:I], Bin)
    @variable(model, t[1:N] >= -1e6)

    @objective(model, Min, data.fixed_costs'* x + fill(1/N, N)' * t)

    max_demand = maximum(sum(demands) for demands in data.demands)
    @constraint(model, capacity, sum(data.capacities[i] * x[i] for i in 1:I) >= max_demand)

    return (x = x, ), t
end

function customize_sub_model!(model::Model, data::SCFLPData, scen_idx::Int; x) 
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, MOI.Silent() => true)

    set_optimizer(model, optimizer)

    I, J = data.n_facilities, data.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    # Set objective
    cost_demands = data.costs .* data.demands[scen_idx]'
    @objective(model, Min, sum(cost_demands .* y))
    # Add constraints
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.demands[scen_idx][:] .* y[i,:]) <= data.capacities[i] * x[i])
end