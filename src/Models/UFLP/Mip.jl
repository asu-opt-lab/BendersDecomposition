struct UFLPMipEnv <: AbstractMipEnv
    model::Model
end

function UFLPMipEnv(data; is_CPLEX=true)

    if is_CPLEX
        model =  Model(CPLEX.Optimizer)
    else
        model =  Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "InfUnbdInfo", 1)
    end
    set_optimizer_attribute(model, MOI.Silent(),true)

    @variable(model, x[1:data.n_facilities], Bin)
    # @variable(model, 0<=x[1:data.n_facilities]<=1)
    @variable(model, y[1:data.n_facilities, 1:data.n_customers] >= 0)

    @objective(model, Min, sum(data.costs[i,j] * data.demands[j] * y[i,j] for i in 1:data.n_facilities, j in 1:data.n_customers) + sum(data.fixed_costs[i] * x[i] for i in 1:data.n_facilities))

    @constraint(model, demand[j in 1:data.n_customers], sum(y[i,j] for i in 1:data.n_facilities) == 1)
    @constraint(model, cons[i in 1:data.n_facilities, j in 1:data.n_customers],   y[i,j] <= x[i])
    
    return UFLPMipEnv(model)
end