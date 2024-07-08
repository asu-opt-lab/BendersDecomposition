
function generate_cut(
    master_env::AbstractMasterEnv, 
    sub_env::AbstractSubEnv,
    ::OrdinaryCutStrategy;
    time_limit::Float64 = 1000.00)

    cVal = master_env.value_x
    for i in eachindex(cVal)
        set_normalized_rhs(sub_env.cconstr[i], cVal[i])
    end
    
    start_time = time()
    set_time_limit_sec(sub_env.model, max(time_limit,60))
    JuMP.optimize!(sub_env.model)
    sub_time = time() - start_time

    status = dual_status(sub_env.model)
    subObjVal = 1e+99

    if status == FEASIBLE_POINT
        subObjVal = JuMP.objective_value(sub_env.model)
        ex = @expression(master_env.model, 
        - master_env.var["t"] + subObjVal 
        + sum(dual(sub_env.cconstr[i]) * (master_env.var["cvar"][i] - cVal[i]) for i in eachindex(cVal)))  

    elseif status == INFEASIBILITY_CERTIFICATE
        if has_duals(sub_env.model)
            subObjVal = JuMP.objective_value(sub_env.model)
            check = sum(dual(sub_env.sub_constr[i]) * sub_env.sub_rhs[i] for i in 1:length(sub_env.sub_rhs)) + sum(dual(sub_env.cconstr[i]) * cVal[i] for i in eachindex(cVal))
            ex = @expression(master_env.model, 
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
        # throw(-1)
        ex = 0
        subObjVal = -Inf
    end         



    sub_env.obj_value = subObjVal
    # @info "subObjVal: $subObjVal"
    return sub_time,ex
end


# function generate_cut(
#     master_env::SAACFLPMasterEnv, 
#     sub_env::AbstractSubEnv,
#     w,
#     ::OrdinaryCutStrategy;
#     time_limit::Float64 = 1000.00)

#     cVal = master_env.value_x
#     for i in eachindex(cVal)
#         set_normalized_rhs(sub_env.cconstr[i], cVal[i])
#     end
    
#     start_time = time()
#     # set_time_limit_sec(sub_env.model, time_limit)
#     JuMP.optimize!(sub_env.model)
#     sub_time = time() - start_time

#     status = dual_status(sub_env.model)
#     subObjVal = 1e+99

#     if status == FEASIBLE_POINT
#         subObjVal = JuMP.objective_value(sub_env.model)
#         ex = @expression(master_env.model, 
#         - master_env.var["t"][w] + subObjVal 
#         + sum(dual(sub_env.cconstr[i]) * (master_env.var["cvar"][i] - cVal[i]) for i in eachindex(cVal)))  

#     elseif status == INFEASIBILITY_CERTIFICATE
#         if has_duals(sub_env.model)
#             subObjVal = JuMP.objective_value(sub_env.model)
#             check = sum(dual(sub_env.sub_constr[i]) * sub_env.sub_rhs[i] for i in 1:length(sub_env.sub_rhs)) + sum(dual(sub_env.cconstr[i]) * cVal[i] for i in eachindex(cVal))
#             ex = @expression(master_env.model, 
#             sum(dual(sub_env.sub_constr[i]) * sub_env.sub_rhs[i] for i in 1:length(sub_env.sub_rhs)) 
#             + sum(dual(sub_env.cconstr[i]) * master_env.var["cvar"][i] for i in 1:length(cVal)))
#             @info "Feasible Cut"
#         else
#             @error "infeasible sub has no infeasible ray"
#             throw(-1)
#         end      
#         subObjVal = 1e+99
#     else
#         @error "dual of sub is neither feasible nor infeasible certificate: $status"
#         throw(-1)
#     end         

#     sub_env.obj_value = subObjVal

#     return sub_time,ex
# end
