export UFLPMasterProblem

mutable struct UFLPMasterEnv <: AbstractMasterEnv
    model::Model
    var::Dict
    coef::Vector{Float64}
    value_x::Vector{Float64}
    value_t::Float64
    obj_value::Float64
    best_solution
    solution_in
end

function UFLPMasterProblem(data::UFLPData; solver::Symbol=:Gurobi)


    if solver == :CPLEX
        model = Model(CPLEX.Optimizer)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        # set_optimizer_attribute(model, "InfUnbdInfo", 1)
    end
    set_optimizer_attribute(model, MOI.Silent(),true)

    # pre
    N = data.n_facilities
    # Variables
    @variable(model, x[1:N], Bin)
    # @variable(model, 0<=x[1:N]<=1)
    @variable(model, t >= -1e06)

    # Objective
    @objective(model, Min, sum(data.fixed_costs[i] * x[i] for i in 1:N) + t)

    return UFLPMasterEnv(model, Dict("cvar"=>x, "t"=>t), data.fixed_costs, zeros(N), 0.0, 0.0, [zeros(N), 0], [ones(N), 1e08])
    
end