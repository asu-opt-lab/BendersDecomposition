export BSPProblem

mutable struct CFLPBSPEnv <: AbstractSubEnv
    model::Model
end

#### Ordinary version
function generate_BSPProblem(data::CFLPData; solver::Symbol=:Gurobi)

    if solver == :CPLEX
        model = Model(CPLEX.Optimizer)
        # set_optimizer_attribute(model, "CPX_PARAM_REDUCE", 0)
        # set_optimizer_attribute(model, "CPX_PARAM_LPMETHOD", CPX_ALG_BARRIER)
        # set_optimizer_attribute(model, "CPX_PARAM_LPMETHOD", CPX_ALG_DUAL)
        # set_optimizer_attribute(model, "CPX_PARAM_BARALG", 1)
        # set_optimizer_attribute(model, "CPX_PARAM_BAROBJRNG", 1e75)
        # set_optimizer_attribute(model, "CPX_PARAM_BAREPCOMP", 1e-6)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        # model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "Method", 2)
        set_optimizer_attribute(model, "InfUnbdInfo", 1)
        # set_optimizer_attribute(model, "LPWarmStart", 0)
    end
    set_optimizer_attribute(model, MOI.Silent(),true)
    # set_time_limit_sec(model, 10)
   
    # pre
    N = data.n_facilities
    M = data.n_customers
    
    # Variables
    @variable(model, 1>=y[1:N,1:M]>=0)
    @variable(model, x[1:N])
    # @variable(model, b)

    # Objective
    obj = @expression(model, sum(data.costs[i,j] * data.demands[j] * y[i,j] for i in 1:N, j in 1:M))
    @objective(model, Min, obj)

    # Constraints
    # @constraint(model, c1[j in 1:M], sum(y[i,j] for i in 1:N) == b)
    # @constraint(model, c2[i in 1:N], sum(data.demands[j] * y[i,j] for j in 1:M) <= data.capacities[i] * x[i])
    # @constraint(model, c3[i in 1:N, j in 1:M], y[i,j] <= x[i])
    # @constraint(model, cx[i in 1:N], x[i] == 0) #x̂[i]
    # @constraint(model, cb, b == 0) #b̂


    @constraint(model, cb[j in 1:M], sum(y[i,j] for i in 1:N) == 1)
    @constraint(model, c2[i in 1:N], sum(data.demands[j] * y[i,j] for j in 1:M) <= data.capacities[i] * x[i])
    @constraint(model, c3[i in 1:N, j in 1:M], y[i,j] <= x[i])
    @constraint(model, cx[i in 1:N], x[i] == 0) #x̂[i]
    # @constraint(model, cb, b == 0) #b̂

    return CFLPBSPEnv(model)
end


