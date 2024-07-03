
function generate_cut(
    master_env::AbstractMasterEnv, 
    sub_env::AbstractSubEnv, 
    ::KNCutStrategy;
    time_limit::Float64 = 1000.00)

    data = sub_env.data
    # knapsack_subproblems = sub_env.knapsack_subproblems
    I = data.n_facilities
    J = data.n_customers

    cVal = master_env.value_x
    for i in eachindex(cVal)
        set_normalized_rhs(sub_env.cconstr[i], cVal[i])
    end

    time_start = time()
    set_time_limit_sec(master_env.model, time_limit)
    optimize!(sub_env.model)
    time_end = time()-time_start
    subObjVal = objective_value(sub_env.model)
    u = dual.(sub_env.model[:c1])
    KP = zeros(I)

    for i in 1:I
        model = sub_env.knapsack_subproblems[i]
        @objective(model, Min, sum((data.demands[j] * data.costs[i,j] - u[j]) * model[:z][j] for j in 1:J))
        optimize!(model)
        KP[i] = objective_value(model)
    end

    ex = @expression(master_env.model, - master_env.var["t"] + sum(u) + sum(KP[i] * master_env.var["cvar"][i] for i in 1:I))

    sub_env.obj_value = subObjVal

    return time_end,ex
end