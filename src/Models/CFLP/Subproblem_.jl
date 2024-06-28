
mutable struct CFLPBSPADEnv <: AbstractSubEnv
    model::Model
    constr::Array
    rhs::Array
    cconstr::Array
    obj::AffExpr
    oconstr::Any
end

function generate_BSPProblem_Advanced(data::CFLPData; solver::Symbol=:CPLEX)

    if solver == :CPLEX
        model = Model(CPLEX.Optimizer)
        # set_optimizer_attribute(model, "CPX_PARAM_REDUCE", 0)
    elseif solver == :Gurobi
        model = Model(Gurobi.Optimizer)
        # set_optimizer_attribute(model, "Method", 1)
        # set_optimizer_attribute(model, "InfUnbdInfo", 1)
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

    @constraint(model, cb[j in 1:M], sum(y[i,j] for i in 1:N) == 1)
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


    return CFLPBSPADEnv(model, constr, rhs, cconstr,obj, oconstr)
end