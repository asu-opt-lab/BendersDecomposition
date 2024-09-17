export UFLPSplitSubEnv

mutable struct UFLPStandardSubEnv <: AbstractSubEnv
    model::Model
    sub_constr::Array
    sub_rhs::Array
    cconstr::Array
    obj_value::Float64

    algo_params
end

mutable struct UFLPStandardADSubEnv <: AbstractSubEnv
    model::Model
    sub_constr::Array
    sub_rhs::Array
    cconstr::Array
    obj_value::Float64
    obj
    oconstr

    algo_params
end

mutable struct UFLPStandardDualADSubEnv <: AbstractSubEnv
    model::Model
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


function UFLPStandardSubEnv(data, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr = generate_UFLP_subproblem(data, solver=solver) 
    # model, constr, rhs, cconstr = generate_UFLP_subproblem_Advanced(data, solver=solver) 

    return UFLPStandardSubEnv(model, constr, rhs, cconstr, 0.0, algo_params)
end

function UFLPStandardADSubEnv(data, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr, obj, oconstr = generate_UFLP_subproblem_Advanced(data, solver=solver) 

    return UFLPStandardADSubEnv(model, constr, rhs, cconstr, 0.0, obj, oconstr, algo_params)

    # model = generate_UFLP_dualsubproblem_Advanced(data, solver=solver) 

    return UFLPStandardDualADSubEnv(model, 0.0,algo_params)
end

function UFLPSplitSubEnv(data, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr = generate_UFLP_subproblem(data; solver=solver) 

    BSPProblem = generate_BSPProblem(data; solver=solver)
    BSPProblem2 = generate_BSPProblem(data; solver=solver)

    # BSPProblem = generate_BSPProblem_Advanced(data; solver=solver)
    # BSPProblem2 = generate_BSPProblem_Advanced(data; solver=solver)

    # BSPProblem = generate_UFLP_dualsubproblem_Advanced(data; solver=solver)
    # BSPProblem2 = generate_UFLP_dualsubproblem_Advanced(data; solver=solver)
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

# function generate_UFLP_subproblem_Advanced(data::UFLPData; solver::Symbol=:Gurobi)
function generate_UFLP_subproblem_Advanced(data; solver::Symbol=:Gurobi)
    if solver == :CPLEX
        model = Model(CPLEX.Optimizer)
        # set_optimizer_attribute(model, "CPX_PARAM_REDUCE", 0)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "Method", 2)
        # set_optimizer_attribute(model, "InfUnbdInfo", 1)
    end
    set_optimizer_attribute(model, MOI.Silent(),true)
    # set_time_limit_sec(model, 10)
    # pre
    N = data.n_facilities
    M = data.n_customers
    
    @variable(model, y[1:data.n_facilities, 1:data.n_customers] >= 0)

    cvar = Dict()
    cvar["x"] = @variable(model, x[1:data.n_facilities])
    @variable(model, σ)

    @objective(model, Min, σ)

    @constraint(model, cb[j in 1:M], sum(y[i,j] for i in 1:N) + σ >= 1)
    @constraint(model, cbb[j in 1:M], -sum(y[i,j] for i in 1:N) + σ >= -1)
    # @constraint(model, c2[i in 1:N], sum(data.demands[j] * y[i,j] for j in 1:M) <= data.capacities[i] * cvar["x"][i])
    @constraint(model, c3[i in 1:N, j in 1:M], -y[i,j] + σ >= -cvar["x"][i])
    # @constraint(model, c3[i in 1:N, j in 1:M], -y[i,j] + σ >= 0)
    constr = []
    rhs = []

    # all_affine_constraints = [all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.GreaterThan{Float64});
    #                           all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64});
    #                           all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.EqualTo{Float64})]

    # for (i, c) in enumerate(all_affine_constraints)
    #     push!(constr, c)
    #     push!(rhs, normalized_rhs(c))
    # end

    cconstr = []

    for i in 1:data.n_facilities
        push!(cconstr, @constraint(model, cvar["x"][i] == 0)) #x̂
    end


    obj = @expression(model, sum(data.costs[i,j] * data.demands[j] * y[i,j] for i in 1:data.n_facilities, j in 1:data.n_customers))
    oconstr = @constraint(model, -obj + σ >= 0) #-η


    return model, constr, rhs, cconstr, obj, oconstr
end





function generate_UFLP_dualsubproblem_Advanced(data; solver::Symbol=:Gurobi)
    if solver == :CPLEX
        model = Model(CPLEX.Optimizer)
        # set_optimizer_attribute(model, "CPX_PARAM_REDUCE", 0)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "Method", 2)
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

    @constraint(model, sum(π1) + sum(π2) + sum(π3) + π0 == 1)

    return model
end


