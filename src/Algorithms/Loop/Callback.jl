function run_Benders_callback(
    data::AbstractData,
    master_env::AbstractMasterEnv,
    sub_env::AbstractSubEnv)
    
    tic = time()
    set_attribute(master_env.model, MOI.LazyConstraintCallback(), my_callback)

    global Master_env = master_env
    global Sub_env = sub_env
    global Data = data
    global number_of_subproblem_solves = 0

    JuMP.optimize!(master_env.model)
    toc = time()
    @info "Time to compute objective value: $(toc - tic)"
    return JuMP.objective_value(master_env.model)
end

function my_callback(cb_data)
    status = JuMP.callback_node_status(cb_data, Master_env.model)
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER
        global number_of_subproblem_solves += 1
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

