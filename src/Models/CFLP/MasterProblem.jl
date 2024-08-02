export MasterProblem

mutable struct CFLPMasterEnv <: AbstractMasterEnv
    model::Model
    var::Dict
    coef::Vector{Int}
    value_x::Vector{Float64}
    value_t::Float64
    obj_value::Float64
end

function CFLPMasterProblem(data::CFLPData; solver::Symbol=:Gurobi)


    if solver == :CPLEX
        model = Model(CPLEX.Optimizer)
        # model = direct_model(CPLEX.Optimizer())
        # set_optimizer_attribute(model, "CPX_PARAM_EPINT", 1e-03)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        # set_optimizer_attribute(model, "InfUnbdInfo", 1)
    end
    set_optimizer_attribute(model, MOI.Silent(),true)
    # set_time_limit_sec(model, 10)

    println("########### Building Master Problem ###########")
    # pre
    N = data.n_facilities
    # Variables
    @variable(model, x[1:N], Bin)
    # @variable(model, 0<=x[1:N]<=1)
    @variable(model, t >= -1e06)

    # Objective
    @objective(model, Min, sum(data.fixed_costs[i] * x[i] for i in 1:N) + t)

    # Constraints
    
    @constraint(model, sum(data.capacities[i] * x[i] for i in 1:N) >= sum(data.demands))
    

    return CFLPMasterEnv(model, Dict("cvar"=>x, "t"=>t), data.fixed_costs, zeros(N), 0.0, 0.0)
    
end
