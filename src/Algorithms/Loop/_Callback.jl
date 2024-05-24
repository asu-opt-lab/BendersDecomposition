function run_Benders_callback(
    data::AbstractData,
    master_env::AbstractMasterEnv,
    sub_env::AbstractSubEnv)
    
    tic = time()
    set_attribute(master_env.model, MOI.LazyConstraintCallback(), lazy_callback)
    # set_attribute(master_env.model, MOI.UserCutCallback(), user_callback)
    # MOI.set(master_env.model, Gurobi.CallbackFunction(), lazy_callback)
    global Master_env = master_env
    global Sub_env = sub_env
    global Data = data
    global number_of_subproblem_solves = 0
    global explored_nodes = []
    global unexplored_nodes = []
    global best_upper_bound = []

    global cb_calls = Cint[]
    JuMP.optimize!(master_env.model)
    toc = time()
    @info "Time to compute objective value: $(toc - tic)"
    
    x_col = Gurobi.c_column(master_env.model, master_env.model[:t])
    GRBgetdblattrelement(master_env.model, "LB", x_col, Ref{Cdouble}(NaN))
    return JuMP.objective_value(master_env.model)
end

function lazy_callback(cb_data)

    status = JuMP.callback_node_status(cb_data, Master_env.model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER
        global number_of_subproblem_solves += 1
        # @info "add"
        # @info "number_of_subproblem_solves = $number_of_subproblem_solves"
        Master_env.value_x = JuMP.callback_value.(cb_data, Master_env.var["cvar"])

        Master_env.value_t = JuMP.callback_value.(cb_data, Master_env.var["t"])
        _,ex = generate_cut(Master_env, Sub_env, Sub_env.algo_params.cut_strategy)
        # callback_generate_cut(Master_env, Subb_env, CVal, Loop_strategy, Cut_strategy, Data, Knapsack_subproblems, cb_data, EtaVal_new)
        cons = @build_constraint(0>=ex)
        MOI.submit(Master_env.model, MOI.LazyConstraint(cb_data), cons)
        
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


