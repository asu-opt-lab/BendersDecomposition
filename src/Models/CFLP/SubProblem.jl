export SubProblem

mutable struct CFLPStandardSubEnv <: AbstractSubEnv
    model::Model
    sub_constr::Array
    sub_rhs::Array
    cconstr::Array
    obj_value::Float64

    algo_params
end

mutable struct CFLPStandardADSubEnv <: AbstractSubEnv
    model::Model
    sub_constr::Array
    sub_rhs::Array
    cconstr::Array
    obj_value::Float64
    obj::AffExpr
    oconstr::Any

    algo_params
end


mutable struct CFLPSplitSubEnv <: AbstractSubEnv
    model::Model
    sub_constr::Array
    sub_rhs::Array
    cconstr::Array
    obj_value::Float64

    algo_params
    data
    BSPProblem
    split_info
    BSPProblem2
end


function CFLPStandardSubEnv(data::CFLPData, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr = generate_CFLP_subproblem(data, solver=solver) 

    return CFLPStandardSubEnv(model, constr, rhs, cconstr, 0.0, algo_params)
end



function CFLPStandardADSubEnv(data::CFLPData, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr, obj, oconstr = generate_CFLP_subproblem_Advanced(data, solver=solver) 

    return CFLPStandardADSubEnv(model, constr, rhs, cconstr, 0.0, obj, oconstr, algo_params)
end

function CFLPSplitSubEnv(data::CFLPData, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr = generate_CFLP_subproblem(data; solver=solver) 

    BSPProblem = generate_BSPProblem(data; solver=solver)
    BSPProblem2 = generate_BSPProblem(data; solver=solver)
   
    split_info = SplitInfo([],[],[])
    
    return CFLPSplitSubEnv(model, constr, rhs, cconstr, 0.0, algo_params, data, BSPProblem, split_info, BSPProblem2)
end




mutable struct CFLPStandardKNSubEnv <: AbstractSubEnv
    model::Model
    sub_constr::Array
    sub_rhs::Array
    cconstr::Array
    obj_value::Float64

    algo_params
    data::CFLPData
    knapsack_subproblems
end

function CFLPStandardKNSubEnv(data::CFLPData, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr = generate_CFLP_subproblem(data, solver=solver) 

    kanpsack_subproblems = generate_knapsack_subproblems(data)
    return CFLPStandardKNSubEnv(model, constr, rhs, cconstr, 0.0, algo_params, data, kanpsack_subproblems)
end

#### Knapsack version
function generate_knapsack_subproblems(data::CFLPData)
    I = data.n_facilities
    J = data.n_customers
    models = Vector{Model}(undef, I)
    u = zeros(J) 

    for i in 1:I
        model = Model(CPLEX.Optimizer)
        set_optimizer_attribute(model, MOI.Silent(),true)
        @variable(model, 0<=z[1:data.n_customers]<=1)
        @objective(model, Min, sum((data.demands[j] * data.costs[i,j] - u[j]) * z[j] for j in 1:J)) 

        @constraint(model, sum(data.demands[j] * z[j] for j in 1:J) <= data.capacities[i])

        models[i] = model
    end
    return models
end