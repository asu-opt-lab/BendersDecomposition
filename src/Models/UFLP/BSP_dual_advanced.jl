mutable struct UFLPBSPADEnv <: AbstractSubEnv
    model::Model
end

function generate_UFLP_dualsubproblem_Advanced(data; solver::Symbol=:Gurobi)
    if solver == :CPLEX
        model = Model(CPLEX.Optimizer)
        # set_optimizer_attribute(model, "CPX_PARAM_REDUCE", 0)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        # set_optimizer_attribute(model, "Method", 2)
        # set_optimizer_attribute(model, "InfUnbdInfo", 1)
    end
    set_optimizer_attribute(model, MOI.Silent(),true)
    # set_time_limit_sec(model, 10)
    # pre
    N = data.n_facilities
    M = data.n_customers
    

    @variable(model, π1[1:M]>=0)
    @variable(model, π2[1:M]>=0)
    @variable(model, π3[1:N,1:M]>=0)
    @variable(model, π0>=0)


    # @objective(model, Min, sum(π1) - sum(π2[i,j] for i in 1:N, j in 1:M) + π0)

    @constraint(model, con[i=1:N, j=1:M], π1[j] - π2[j] - π3[i,j] <= data.costs[i,j] * data.demands[j] * π0)

    # @constraint(model, sum(π1) + sum(π2) + sum(π3) + π0 == 1)

    @constraint(model, sum(π3) + π0 == 1)

    return UFLPBSPADEnv(model)
end