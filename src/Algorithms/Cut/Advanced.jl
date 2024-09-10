

# function generate_cut(
#     master_env::AbstractMasterEnv, 
#     sub_env::AbstractSubEnv, 
#     ::AdvancedCutStrategy;
#     time_limit::Float64 = 1000.00)

#     cVal = master_env.value_x
#     l = length(cVal)    
#     for i in eachindex(cVal)
#         set_normalized_rhs(sub_env.cconstr[i], cVal[i])
#         set_normalized_rhs(sub_env.cconstr[l+i], -cVal[i])
#     end
#     set_normalized_rhs(sub_env.oconstr, -master_env.value_t)

#     # set_time_limit_sec(sub_env.model, max(time_limit,100))
#     start_time = time()
#     JuMP.optimize!(sub_env.model)
#     sub_time = time() - start_time

#     status = dual_status(sub_env.model)
#     # update value_η!!!!!!!!!!!!!!!!!!!
#     if status == FEASIBLE_POINT
#         subObjVal = JuMP.objective_value(sub_env.model) 
#         ex = @expression(master_env.model,  
#         subObjVal 
#         - dual(sub_env.oconstr)*(master_env.var["t"] - master_env.value_t)
#         - sum(dual(sub_env.cconstr[i]) * (master_env.var["cvar"][i] - cVal[i]) for i in eachindex(cVal)) 
#         - sum(dual(sub_env.cconstr[l+i]) * (master_env.var["cvar"][i] - cVal[i]) for i in eachindex(cVal)))       
#     else
#         @error "dual of sub is neither feasible nor infeasible certificate: $status"
#         throw(-1)
#     end

#     if abs(subObjVal) <= 1e-06
#         sub_env.obj_value = value(sub_env.obj)#sum(data.costs[i,j] * data.demands[j] * JuMP.value(sub_env.model[:y][i,j]) for i in 1:data.n_facilities, j in 1:data.n_customers)  
#     else 
#         sub_env.obj_value = Inf
#     end

#     return sub_time,ex
# end




function generate_cut(
    master_env::AbstractMasterEnv, 
    sub_env::AbstractSubEnv, 
    ::AdvancedCutStrategy;
    time_limit::Float64 = 1000.00)

    cVal = master_env.value_x
    l = length(cVal)    
    # for i in eachindex(cVal)
    #     set_normalized_rhs(sub_env.cconstr[i], cVal[i])
    # end
    for i in length(cVal)
        for j in length(sub_env.model[:cb])
            set_normalized_rhs(sub_env.model[:c3][i,j], -cVal[i])
        end
    end
    # set_normalized_rhs.(sub_env.model[:c3], -cVal)
    set_normalized_rhs(sub_env.oconstr, -master_env.value_t)

    # set_time_limit_sec(sub_env.model, max(time_limit,100))
    start_time = time()
    JuMP.optimize!(sub_env.model)
    sub_time = time() - start_time

    status = dual_status(sub_env.model)
    # update value_η!!!!!!!!!!!!!!!!!!!
    if status == FEASIBLE_POINT || status == INFEASIBILITY_CERTIFICATE
        subObjVal = JuMP.objective_value(sub_env.model) 
        @info subObjVal
        # ex = @expression(master_env.model,  
        # subObjVal 
        # - dual(sub_env.oconstr)*(master_env.var["t"] - master_env.value_t)
        # + sum(dual(sub_env.cconstr[i]) * (master_env.var["cvar"][i] - cVal[i]) for i in eachindex(cVal)))    
        ex = @expression(master_env.model, -sum(dual.(sub_env.model[:cb])) - sum(dual.(sub_env.model[:c3][i])*master_env.model[:x][i] for i in eachindex(cVal)) - dual(sub_env.oconstr)*master_env.model[:t])
        @info ex
        # @info ex1
    else
        @error "dual of sub is neither feasible nor infeasible certificate: $status"
        throw(-1)
    end

    # if abs(subObjVal) <= 1e-06
    #     sub_env.obj_value = value(sub_env.obj)#sum(data.costs[i,j] * data.demands[j] * JuMP.value(sub_env.model[:y][i,j]) for i in 1:data.n_facilities, j in 1:data.n_customers)  
    # else 
    #     sub_env.obj_value = Inf
    # end
    sub_env.obj_value = Inf
    return sub_time,ex
end




