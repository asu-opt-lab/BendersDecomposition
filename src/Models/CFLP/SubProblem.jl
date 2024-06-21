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
end


function CFLPStandardSubEnv(data::CFLPData, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr = generate_CFLP_subproblem(data, solver=solver) 

    return CFLPStandardSubEnv(model, constr, rhs, cconstr, 0.0, algo_params)
end

function CFLPSplitSubEnv(data::CFLPData, algo_params; solver::Symbol=:Gurobi)
    
    model, constr, rhs, cconstr = generate_CFLP_subproblem(data; solver=solver) 
    BSPProblem = generate_BSPProblem(data; solver=solver)
    split_info = SplitInfo([],[],[])
    
    return CFLPSplitSubEnv(model, constr, rhs, cconstr, 0.0, algo_params, data, BSPProblem, split_info)
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


function generate_CFLP_subproblem(data::CFLPData; solver::Symbol=:CPLEX)
    if solver == :CPLEX
        model =  Model(CPLEX.Optimizer)
        set_optimizer_attribute(model, "CPX_PARAM_REDUCE", 0)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "Method", 1)
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
    @constraint(model, c2[i in 1:N], sum(data.demands[j] * y[i,j] for j in 1:M) <= data.capacities[i] * cvar["x"][i])
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



export BSPProblem

mutable struct CFLPBSPEnv <: AbstractSubEnv
    model::Model
end

function generate_BSPProblem(data::CFLPData; solver::Symbol=:CPLEX)

    if solver == :CPLEX
        model = Model(CPLEX.Optimizer)
        # set_optimizer_attribute(model, "CPX_PARAM_REDUCE", 0)
        # set_optimizer_attribute(model, "CPX_PARAM_LPMETHOD", CPX_ALG_DUAL)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        # model = Model(Ipopt.Optimizer)
        set_optimizer_attribute(model, "Method", 1)
        set_optimizer_attribute(model, "InfUnbdInfo", 1)
        set_optimizer_attribute(model, "LPWarmStart", 0)
    end
    set_optimizer_attribute(model, MOI.Silent(),true)
   
    # pre
    N = data.n_facilities
    M = data.n_customers
    
    # Variables
    @variable(model, y[1:N,1:M]>=0)
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



function generate_CFLP_subproblem_Advanced(data::CFLPData; solver::Symbol=:CPLEX)

    if solver == :CPLEX
        model = Model(CPLEX.Optimizer)
        # set_optimizer_attribute(model, "CPX_PARAM_REDUCE", 0)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "Method", 1)
        set_optimizer_attribute(model, "InfUnbdInfo", 1)
    end
    set_optimizer_attribute(model, MOI.Silent(),true)

    # pre
    N = data.n_facilities
    M = data.n_customers
    
    @variable(model, y[1:data.n_facilities, 1:data.n_customers] >= 0)

    cvar = Dict()
    cvar["x"] = @variable(model, x[1:data.n_facilities])
    @variable(model, σ)

    @objective(model, Min, σ)

    @constraint(model, c1[j in 1:M], sum(y[i,j] for i in 1:N) == 1)
    @constraint(model, c2[i in 1:N], sum(data.demands[j] * y[i,j] for j in 1:M) <= data.capacities[i] * cvar["x"][i])
    @constraint(model, c3[i in 1:N, j in 1:M], y[i,j] <= cvar["x"][i])

    constr = []
    rhs = []

    all_affine_constraints = [all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.GreaterThan{Float64});
                              all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.LessThan{Float64});
                              all_constraints(model, GenericAffExpr{Float64,VariableRef}, MOI.EqualTo{Float64})]

    for (i, c) in enumerate(all_affine_constraints)
        push!(constr, c)
        push!(rhs, normalized_rhs(c))
    end

    cconstr = []

    for i in 1:data.n_facilities
        push!(cconstr, @constraint(model, cvar["x"][i] + σ >= 0)) #x̂
    end

    for i in 1:data.n_facilities
        push!(cconstr, @constraint(model, -cvar["x"][i] + σ >= 0))
    end


    obj = @expression(model, sum(data.costs[i,j] * data.demands[j] * y[i,j] for i in 1:data.n_facilities, j in 1:data.n_customers))
    oconstr = @constraint(model, -obj + σ >= 0) #-η


    return model, constr, rhs, cconstr,obj, oconstr
end