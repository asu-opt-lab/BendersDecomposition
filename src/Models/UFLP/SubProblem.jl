export UFLPSplitSubEnv

mutable struct UFLPStandardSubEnv <: AbstractSubEnv
    model::Model
    sub_constr::Array
    sub_rhs::Array
    cconstr::Array
    obj_value::Float64

    algo_params
end

mutable struct UFLPSplitSubEnv <: AbstractSubEnv
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


function UFLPStandardSubEnv(data::UFLPData, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr = generate_UFLP_subproblem(data, solver=solver) 

    return UFLPStandardSubEnv(model, constr, rhs, cconstr, 0.0, algo_params)
end

function UFLPSplitSubEnv(data::UFLPData, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr = generate_UFLP_subproblem(data; solver=solver) 

    # BSPProblem = generate_BSPProblem(data; solver=solver)
    # BSPProblem2 = generate_BSPProblem(data; solver=solver)

    BSPProblem = generate_BSPProblem_Advanced(data; solver=solver)
    BSPProblem2 = generate_BSPProblem_Advanced(data; solver=solver)
    println("########### Building two Split SubProblems ###########")

    split_info = SplitInfo([],[],[])
    
    return UFLPSplitSubEnv(model, constr, rhs, cconstr, 0.0, algo_params, data, BSPProblem, split_info, BSPProblem2)
end


function generate_UFLP_subproblem(data::Union{UFLPData,CFLPData}; solver::Symbol=:Gurobi)
    if solver == :CPLEX
        model =  Model(CPLEX.Optimizer)
        # set_optimizer_attribute(model, "CPX_PARAM_REDUCE", 0)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        # set_optimizer_attribute(model, "Method", 2)
        set_optimizer_attribute(model, "InfUnbdInfo", 1)
    end
    set_optimizer_attribute(model, MOI.Silent(),true)

    # pre
    N = data.n_facilities
    M = data.n_customers
    
    # Variables
    @variable(model, y[1:N,1:M]>=0)
    # @variable(model, x[1:N])
    cvar = Dict()
    cvar["x"] = @variable(model, x[1:data.n_facilities])

    # Objective
    obj = @expression(model, sum(data.costs[i,j] * data.demands[j] * y[i,j] for i in 1:N, j in 1:M))
    @objective(model, Min, obj)

    # Constraints
    @constraint(model, c1[j in 1:M], sum(y[i,j] for i in 1:N) == 1)
    # @constraint(model, c2[i in 1:N], sum(data.demands[j] * y[i,j] for j in 1:M) <= data.capacities[i] * cvar["x"][i])
    @constraint(model, c3[i in 1:N, j in 1:M], y[i,j] <= cvar["x"][i])

    constr = []
    rhs = []

    all_affine_constraints = [all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.GreaterThan{Float64});
                              all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64});
                              all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.EqualTo{Float64});]

    for (i, c) in enumerate(all_affine_constraints)
        push!(constr, c)
        push!(rhs, normalized_rhs(c))
    end

    cconstr = []

    for i in 1:data.n_facilities
        push!(cconstr, @constraint(model, cvar["x"][i] == 0))
    end

    return model, constr, rhs, cconstr
end






