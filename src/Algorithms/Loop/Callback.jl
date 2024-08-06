function run_Benders_callback(
    data::AbstractData,
    master_env::AbstractMasterEnv,
    sub_env::AbstractSubEnv)
    
    tic = time()
    set_attribute(master_env.model, MOI.LazyConstraintCallback(), lazy_callback)
    # set_attribute(master_env.model, MOI.LazyConstraintCallback(), my_callback_function)
    # set_attribute(master_env.model, MOI.UserCutCallback(), user_callback)
    # set_attribute(master_env.model, MOI.NumberOfThreads(), 1)
    set_time_limit_sec(master_env.model, 300)
    set_optimizer_attribute(master_env.model, MOI.Silent(),false)

    global Master_env = master_env
    global Sub_env = sub_env
    global Data = data
    global number_of_subproblem_solves = 0
    global number_of_splitproblem_solves = 0


    JuMP.optimize!(master_env.model)
    @info termination_status(master_env.model)
    # @info explored_nodes
    # @info unexplored_nodes
    # @info best_upper_bound
    toc = time()

    @info JuMP.node_count(master_env.model)
    @info JuMP.objective_bound(master_env.model)
    @info JuMP.objective_value(master_env.model)
    @info JuMP.relative_gap(master_env.model)

    @info "Time to compute objective value: $(toc - tic)"


    return JuMP.objective_value(master_env.model)
end

function lazy_callback(cb_data)
    status = JuMP.callback_node_status(cb_data, Master_env.model)
    # @info status
    if status == MOI.CALLBACK_NODE_STATUS_INTEGER
        global number_of_subproblem_solves += 1
        global number_of_splitproblem_solves = 0
        Master_env.value_x = JuMP.callback_value.(cb_data, Master_env.var["cvar"])
        Master_env.value_t = JuMP.callback_value.(cb_data, Master_env.var["t"])
        _,ex = generate_cut(Master_env, Sub_env, ORDINARY_CUTSTRATEGY)
        if Master_env.value_t <= Sub_env.obj_value - 1e-06
            cons = @build_constraint(0>=ex)
            MOI.submit(Master_env.model, MOI.LazyConstraint(cb_data), cons)
        end
        
    elseif status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
        @warn "cb status = CALLBACK_NODE_STATUS_UNKNOWN"
    end
end


function user_callback(cb_data)
    status = JuMP.callback_node_status(cb_data, Master_env.model)
    
    depth = Ref{CPXLONG}()
    ret = CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEDEPTH, depth)
    # @info "depth = $(depth[])"
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL && number_of_splitproblem_solves <= 1 && depth[] == 2

        global number_of_splitproblem_solves += 1
        @info "number_of_subproblem_solves = $number_of_subproblem_solves"
        
        Master_env.value_x = JuMP.callback_value.(cb_data, Master_env.var["cvar"])
        Master_env.value_t = JuMP.callback_value.(cb_data, Master_env.var["t"])
        # @info abs.(Master_env.value_x)
        # lb, ub = fill(NaN, Data.n_facilities+1), fill(NaN, Data.n_facilities+1)
        # @assert CPXcallbackgetlocallb(cb_data, lb, 0, length(lb)-1 ) == 0
        # @assert CPXcallbackgetlocalub(cb_data, ub, 0, length(ub)-1 ) == 0
        # @info lb
        # @info ub
        # @info abs.(lb[2:end])
        # @info ub[2:end]
        # println("There are $(count(lb .≈ ub)) fixed variables")
        # a = Int.(lb .≈ ub)[2:end]
        # push!(a, 0)
        # b = Int(count(lb .≈ ub))
        a,b = select_split_set(Master_env, Sub_env.algo_params.SplitSetSelectionPolicy)
        DCGLP_env = DCGLP(Sub_env, a, b, Sub_env.algo_params.SplitCGLPNormType)
        x̂,t̂ = Master_env.value_x, Master_env.value_t
        solve_DCGLP(Master_env,x̂,t̂, DCGLP_env, Sub_env.BSPProblem, Sub_env.BSPProblem2, Sub_env.algo_params.SplitCGLPNormType)
        γ₀, γₓ, γₜ = generate_cut(Master_env, DCGLP_env, Sub_env.algo_params.StrengthenCutStrategy)
        ex = @build_constraint(-γ₀ - γₓ'Master_env.model[:x] - γₜ*Master_env.model[:t] >= 0)
        MOI.submit(Master_env.model, MOI.UserCut(cb_data), ex)
    elseif status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
        @warn "cb status = CALLBACK_NODE_STATUS_UNKNOWN"
    end
end

# function user_callback(cb_data, cb_where)
#     status = JuMP.callback_node_status(cb_data, Master_env.model)
    
#     # depth = Ref{CPXLONG}()
#     # ret = CPXcallbackgetinfolong(cb_data, CPXCALLBACKINFO_NODEDEPTH, depth)
#     # @info "depth = $(depth[])"
#     if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL && number_of_splitproblem_solves <= 1
#         @info "user"
        
#         global number_of_splitproblem_solves += 1
#         @info "number_of_subproblem_solves = $number_of_subproblem_solves"
        
#         Master_env.value_x = JuMP.callback_value.(cb_data, Master_env.var["cvar"])
#         Master_env.value_t = JuMP.callback_value.(cb_data, Master_env.var["t"])
        
#         lb, ub = fill(NaN, Data.n_facilities), fill(NaN, Data.n_facilities)
#         resultP = Ref{Cint}()
#         @info GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_STATUS, resultP)
#         # @assert MOI.get(grb, Gurobi.VariableAttribute("LB"), index.(Master_env.model[:x][1]))
#         # @assert MOI.get(grb, Gurobi.VariableAttribute("LB"), index.(Master_env.model[:x]))
#         @info lb,ub
#         println("There are $(count(lb .≈ ub)) fixed variables")
#         a = Int.(lb .≈ ub)
#         b = Int(count(lb .≈ ub))

    
#         γ₀, γₓ, γₜ = generate_cut_callback(Master_env, Sub_env, a, b, SPLIT_CUTSTRATEGY)
#         ex = @build_constraint(-γ₀ - γₓ'Master_env.model[:x] - γₜ*Master_env.model[:t] >= 0)
#         MOI.submit(Master_env.model, MOI.UserCut(cb_data), ex)
#         @info "user end"
#     elseif status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
#         @warn "cb status = CALLBACK_NODE_STATUS_UNKNOWN"
#     end
# end

function generate_cut_callback(
    master_env::AbstractMasterEnv,
    sub_env::CFLPSplitSubEnv,
    a_fixed,
    b_fixed,
    ::SplitCutStrategy;
    time_limit = 1000.00)

    # select the index 
    a,b = select_split_set(master_env, sub_env.algo_params.SplitSetSelectionPolicy)
    # a .+= a_fixed
    # b += b_fixed
    DCGLP_env = DCGLP(sub_env, a, b, sub_env.algo_params.SplitCGLPNormType)

    start_time = time()
    solve_DCGLP(master_env, DCGLP_env, sub_env.BSPProblem, sub_env.algo_params.SplitCGLPNormType; time_limit)
    DCGLP_time = time() - start_time
    
    # add split cut
    γ₀, γₓ, γₜ = generate_cut(master_env, DCGLP_env, sub_env.algo_params.StrengthenCutStrategy)
    # @constraint(master_env.model, -γ₀ - γₓ'master_env.model[:x] - γₜ*master_env.model[:t] >= 0) 
    # push!(sub_env.split_info.γ₀s, γ₀)
    # push!(sub_env.split_info.γₓs, γₓ)
    # push!(sub_env.split_info.γₜs, γₜ)

    # # add benders cut
    # if DCGLP_env.ifsolved == true
    #     generate_cut!(master_env, DCGLP_env, sub_env.algo_params.SplitBendersStrategy)
    # else
    #     generate_cut!(master_env, DCGLP_env, ALL_SPLIT_BENDERS_STRATEGY)
    # end
    
    # sub_time,ex = generate_cut(master_env, sub_env, ORDINARY_CUTSTRATEGY)
    
    return γ₀, γₓ, γₜ
end
