export customize_master_model!, customize_sub_model!, customize_mip_model!

function customize_mip_model!(model::Model, problem::UFLPData)
    
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)

    set_optimizer(model, optimizer)
    
    @variable(model, x[1:problem.n_facilities], Bin)
        
    I, J = problem.n_facilities, problem.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    @variable(model, t[1:J] >= 0)
    
    cost_demands = problem.costs .* problem.demands'
    @objective(model, Min, problem.fixed_costs'* x + sum(t))
    
    @constraint(model, obj[j in 1:J], t[j] >= sum(cost_demands[:,j] .* y[:,j]))
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
end

function customize_master_model!(model::Model, problem::UFLPData)

    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPINT" => 1e-9, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPGAP" => 1e-6, MOI.Silent() => true)

    set_optimizer(model, optimizer)

    @variable(model, x[1:problem.n_facilities], Bin)
    @variable(model, t >= -1e6)
    
    @constraint(model, sum(x) >= 2)
    
    @objective(model, Min, problem.fixed_costs'* x + 1.0 * t)
    
    return (x = x, ), t
end

function customize_sub_model!(model::Model, problem::UFLPData; x) 

    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, "CPX_PARAM_EPRHS" => 1e-9, "CPX_PARAM_EPOPT" => 1e-9, "CPX_PARAM_NUMERICALEMPHASIS" => 1, MOI.Silent() => true)

    set_optimizer(model, optimizer)

    I, J = problem.n_facilities, problem.n_customers
    
    @variable(model, y[1:I, 1:J] >= 0)

    cost_demands = problem.costs .* problem.demands'
    @objective(model, Min, sum(cost_demands .* y))

    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open[j in 1:J], y[:, j] .<= x)
end

# """
# To-do
# """
# function update_model!(oracle::DisjunctiveOracle, data::Data{<:UFLPData})
#     dcglp = oracle.dcglp 

#     @constraint(dcglp, [i=1:2], sum(dcglp[:omega_x][i,:]) >= 2 * dcglp[:omega_0][i])
# end