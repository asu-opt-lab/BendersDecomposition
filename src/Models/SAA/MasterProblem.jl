export MasterProblem

mutable struct SAACFLPMasterEnv <: AbstractMasterEnv
    model::Model
    var::Dict
    coef::Vector{Int}
    value_x::Vector{Float64}
    value_t::Any
    obj_value::Float64
end

function MasterProblem(datas::Array{CFLPData}, num_scenario::Int; solver::Symbol=:Gurobi)


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
    # pre
    N = datas[1].n_facilities
    # Variables
    @variable(model, x[1:N], Bin)
    # @variable(model, 0<=x[1:N]<=1)
    @variable(model, t[1:num_scenario] >= -1e06)

    # Objective
    @objective(model, Min, sum(datas[1].fixed_costs[i] * x[i] for i in 1:N) + mean(t))

    # Constraints
    
    @constraint(model, [w in 1:num_scenario], sum(datas[1].capacities[i] * x[i] for i in 1:N) >= sum(datas[w].demands))
    

    return SAACFLPMasterEnv(model, Dict("cvar"=>x, "t"=>t), datas[1].fixed_costs, zeros(N), zeros(num_scenario), 0.0)
    
end

