export SubProblem

# mutable struct CFLPSubEnv <: AbstractSubEnv
#     model::Model
#     sub_constr::Array
#     sub_rhs::Array
#     cconstr::Array
#     obj_value::Float64
# end

mutable struct CFLPStandardSubEnv <: AbstractSubEnv
    model::Model
    sub_constr::Array
    sub_rhs::Array
    cconstr::Array
    obj_value::Float64

    algo_params
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

mutable struct SplitInfo 
    # indices
    γ₀s
    γₓs
    γₜs
    # ifaddall::Bool
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

function CFLPSplitSubEnv(data::CFLPData, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr = generate_CFLP_subproblem(data; solver=solver) 
    BSPProblem = generate_BSPProblem(data; solver=solver)
    BSPProblem2 = generate_BSPProblem(data; solver=solver)
    split_info = SplitInfo([],[],[])
    

    # BSPProblem = generate_BSPProblem_Advanced(data; solver=solver)
    # BSPProblem2 = generate_BSPProblem_Advanced(data; solver=solver)
    return CFLPSplitSubEnv(model, constr, rhs, cconstr, 0.0, algo_params, data, BSPProblem, split_info, BSPProblem2)
end

function CFLPStandardKNSubEnv(data::CFLPData, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr = generate_CFLP_subproblem(data, solver=solver) 

    kanpsack_subproblems = generate_knapsack_subproblems(data)
    return CFLPStandardKNSubEnv(model, constr, rhs, cconstr, 0.0, algo_params, data, kanpsack_subproblems)
end

function CFLPStandardADSubEnv(data::CFLPData, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr, obj, oconstr = generate_CFLP_subproblem_Advanced(data, solver=solver) 

    return CFLPStandardADSubEnv(model, constr, rhs, cconstr, 0.0, obj, oconstr, algo_params)
end








