export customize_master_model!, customize_sub_model!, customize_mip_model!

function customize_mip_model!(model::Model, problem::CFLPData)

    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)

    set_optimizer(model, optimizer)
    
    I, J = problem.n_facilities, problem.n_customers
    @variable(model, x[1:I], Bin)
    @variable(model, y[1:I, 1:J] >= 0)
    @variable(model, t)
    
    cost_demands = problem.costs .* problem.demands'
    @objective(model, Min, problem.fixed_costs'* x + t)

    # Add constraints
    @constraint(model, t >= sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(problem.demands .* y[i,:]) <= problem.capacities[i] * x[i])
    @constraint(model, capacity_total, sum(problem.capacities[i] * x[i] for i in 1:I) >= sum(problem.demands))
end

function customize_master_model!(model::Model, problem::CFLPData)
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)

    set_optimizer(model, optimizer)
    
    @variable(model, x[1:problem.n_facilities], Bin)
    @variable(model, t >= -1e6)

    @objective(model, Min, problem.fixed_costs'* x + t)

    I = problem.n_facilities
    @constraint(model, capacity, sum(problem.capacities[i] * x[i] for i in 1:I) >= sum(problem.demands))

    return (x = x, ), t
end

function customize_sub_model!(model::Model, problem::CFLPData, scen_idx::Int; x) 
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, MOI.Silent() => true)

    set_optimizer(model, optimizer)

    I, J = problem.n_facilities, problem.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    # Set objective
    cost_demands = problem.costs .* problem.demands'
    @objective(model, Min, sum(cost_demands .* y))
    # Add constraints
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(problem.demands[:] .* y[i,:]) <= problem.capacities[i] * x[i])
end