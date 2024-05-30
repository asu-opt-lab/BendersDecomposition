function run_Benders(
    data::AbstractData,
    master_env::AbstractMasterEnv,
    sub_env::AbstractSubEnv,
    time_limit = 3600)
    
    # Initialize
    UB = Inf
    LB = -Inf
    iter = 1

    algo_start_time = time()
    remaining_time = time_limit

    df = DataFrame(iter = Int[], LB = Float64[], UB = Float64[], gap = Float64[], master_time = Float64[], sub_time = Float64[])

    while true

        #### Master Part ####
        master_time = solve_master!(master_env; time_limit = remaining_time)
        LB = master_env.obj_value
        
        #### Sub Part ####
        remaining_time -= master_time
        sub_time,ex = generate_cut(master_env, sub_env, sub_env.algo_params.cut_strategy; time_limit = remaining_time)

        # Update Parameters
        UB_temp = sum(master_env.coef[i] * master_env.value_x[i] for i in eachindex(master_env.value_x))
        UB_temp += sub_env.obj_value    
        UB = min(UB, UB_temp)
        Gap = 100 * (UB - LB)/ abs(UB) 

        # Store Data
        new_row = (iter, LB, UB, Gap, master_time, sub_time)
        push!(df, new_row)  
        
        # Add Cut
        @constraint(master_env.model, 0 >= ex)   

        # Print
        @printf "%5d     %10.2f   %10.2f    %10.2f  %10.2f  %10.2f\n" iter LB UB Gap master_time sub_time  

        # Stopping Criteria
        # Gap 
        if Gap < 1e-3
            break
        end 

        # Time Limit
        algo_run_time = time()
        spend_time = algo_run_time - algo_start_time
        remaining_time = time_limit - spend_time         
        if spend_time > time_limit 
            @info "Time limit $time_limit reached"
            master_time = solve_master!(master_env; time_limit = remaining_time)
            LB = master_env.obj_value
            new_row = (iter+1, LB, Inf, Inf, master_time, Inf)
            push!(df, new_row)  
            break
        end 

        iter += 1
    end

   
    return df
end



function solve_master!(master_env::AbstractMasterEnv; time_limit=1000)

    start_time = time()
    # set_time_limit_sec(master_env.model, time_limit)
    JuMP.optimize!(master_env.model)
    master_time = time() - start_time

    master_env.obj_value = JuMP.objective_value(master_env.model)
    master_env.value_t = JuMP.value.(master_env.var["t"])
    master_env.value_x = JuMP.value.(master_env.var["cvar"])

    return master_time
end