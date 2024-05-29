function run_Benders_callback(
    data::AbstractData,
    master_env::AbstractMasterEnv,
    sub_env::AbstractSubEnv)
    
    tic = time()
    set_attribute(master_env.model, MOI.LazyConstraintCallback(), lazy_callback)
    # set_attribute(master_env.model, MOI.UserCutCallback(), user_callback)

    global Master_env = master_env
    global Sub_env = sub_env
    global Data = data
    global number_of_subproblem_solves = 0
    global explored_nodes = []
    global unexplored_nodes = []
    global best_upper_bound = []


    JuMP.optimize!(master_env.model)
    @info termination_status(master_env.model)
    # @info explored_nodes
    # @info unexplored_nodes
    # @info best_upper_bound
    toc = time()
    # cpx = JuMP.unsafe_backend(master_env.model)
    # gap = CPXgetmiprelgap(cpx.env, cpx.lp, Ref{Cdouble}())
    # obj = CPXgetobjval(cpx.env, cpx.lp, Ref{Cdouble}())
    @info JuMP.node_count(master_env.model)
    @info JuMP.objective_bound(master_env.model)
    @info JuMP.objective_value(master_env.model)
    @info JuMP.relative_gap(master_env.model)
    # println("Gap: ", gap[])
    # println("Objective value: ", obj[])
    # @info CPXgetnumusercuts(cpx.env, cpx.lp)
    @info "Time to compute objective value: $(toc - tic)"


    return JuMP.objective_value(master_env.model)
end

function lazy_callback(cb_data)
    status = JuMP.callback_node_status(cb_data, Master_env.model)
    # valueP = Ref{Cdouble}()
    # ret = CPXcallbackgetinfodbl(cb_data, CPXCALLBACKINFO_BEST_BND, valueP)
    # @info "Best bound is currently: $(valueP[])"
    # n1 = Ref{CPXLONG}()
    # CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODECOUNT, n1)
    # n2 = Ref{Cdouble}()
    # CPXcallbackgetinfodbl(cb_data, CPXCALLBACKINFO_BEST_BND, n2)
    # n3 = Ref{CPXLONG}()
    # CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODESLEFT, n3)
    # push!(explored_nodes, n1[])
    # push!(best_upper_bound, n2[])
    #     # @info valueP[]
    # push!(unexplored_nodes, n3[])
    
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER
        global number_of_subproblem_solves += 1
        # @info "add"
        
        # @info "number_of_subproblem_solves = $number_of_subproblem_solves"
        Master_env.value_x = JuMP.callback_value.(cb_data, Master_env.var["cvar"])

        Master_env.value_t = JuMP.callback_value.(cb_data, Master_env.var["t"])

        _,ex = generate_cut(Master_env, Sub_env, Sub_env.algo_params.cut_strategy)
        if Master_env.value_t <= Sub_env.obj_value - 1e-06
            cons = @build_constraint(0>=ex)
            MOI.submit(Master_env.model, MOI.LazyConstraint(cb_data), cons)
        end
        # callback_generate_cut(Master_env, Sub_env, Master_env.value_x,  Data, cb_data, Master_env.value_t)
        
    elseif status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
        @warn "cb status = CALLBACK_NODE_STATUS_UNKNOWN"
    end
end

function user_callback(cb_data)
    status = JuMP.callback_node_status(cb_data, Master_env.model)
    # @info "user"
    n = Ref{CPXLONG}()
    CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODECOUNT, n)
    # println(n[])
    
    depth = Ref{CPXLONG}()
    ret = CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEDEPTH, depth)
    # println(depth[])
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL && depth[] == 0 && number_of_splitproblem_solves <= 10
        # @info "user"
        push!(explored_nodes, n[])
        global number_of_splitproblem_solves += 1
        # @info "number_of_subproblem_solves = $number_of_subproblem_solves"
        
        CVal = JuMP.callback_value.(cb_data, Master_env.var["cvar"])
        EtaVal_new = JuMP.callback_value.(cb_data, Master_env.var["t"])
        
        γ₀, γₓ, γₜ = generate_user_cut(Master_env, BSP, CVal, SPLIT_CUTSTRATEGY,Data)
        ex = @build_constraint(-γ₀ - γₓ'Master_env.model[:x] - γₜ*Master_env.model[:t] >= 0)
        MOI.submit(Master_env.model, MOI.UserCut(cb_data), ex)
        @info "user end"
    elseif status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
        @warn "cb status = CALLBACK_NODE_STATUS_UNKNOWN"
    end
end



function callback_generate_cut(
    master_env::AbstractMasterEnv, 
    sub_env::AbstractSubEnv, 
    cVal::Vector{Float64}, 
    data::CFLPData,
    cb_data::Any,
    tval::Any)

    I = data.n_facilities
    J = data.n_customers

    for i in eachindex(cVal)
        set_normalized_rhs(sub_env.cconstr[i], cVal[i])
    end

    optimize!(sub_env.model)
    subObjVal = objective_value(sub_env.model)


    if tval >= subObjVal - 1e-06
        return 
    end

    for i in eachindex(cVal)
        set_normalized_rhs(sub_env.cconstr[i], cVal[i])
    end

    JuMP.optimize!(sub_env.model)
    status = dual_status(sub_env.model)
    subObjVal = 1e+99

    if status == FEASIBLE_POINT
        subObjVal = JuMP.objective_value(sub_env.model)
        ex = @build_constraint( 
        master_env.var["t"] >= subObjVal 
        + sum(dual(sub_env.cconstr[i]) * (master_env.var["cvar"][i] - cVal[i]) for i in eachindex(cVal)))  
        # @info "Optimal Cut"
    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(sub_env.model)
            subObjVal = JuMP.objective_value(sub_env.model)
            check = sum(dual(sub_env.sub_constr[i]) * sub_env.sub_rhs[i] for i in 1:length(sub_env.sub_rhs)) + sum(dual(sub_env.cconstr[i]) * cVal[i] for i in eachindex(cVal))
            ex = @build_constraint(0>= 
            sum(dual(sub_env.sub_constr[i]) * sub_env.sub_rhs[i] for i in 1:length(sub_env.sub_rhs)) 
            + sum(dual(sub_env.cconstr[i]) * master_env.var["cvar"][i] for i in 1:length(cVal)))
            @info "Feasible Cut"
        else
            @error "infeasible sub has no infeasible ray"
            throw(-1)
        end      
        subObjVal = 1e+99
    else
        @error "dual of sub is neither feasible nor infeasible certificate: $status"
        throw(-1)
    end       

    MOI.submit(master_env.model, MOI.LazyConstraint(cb_data), ex)



    return subObjVal
end