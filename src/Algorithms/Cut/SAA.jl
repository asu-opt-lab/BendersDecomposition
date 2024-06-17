function generate_cut(
    master_env::AbstractMasterEnv,
    sub_env::CFLPSplitSubEnv,
    w,
    ::SplitCutStrategy;
    time_limit = 1000.00)

    # select the index 
    a,b = select_split_set(master_env, sub_env.algo_params.SplitSetSelectionPolicy)

    DCGLP_env = DCGLP(sub_env, a, b, sub_env.algo_params.SplitCGLPNormType)

    x̂,t̂ = master_env.value_x, master_env.value_t[w]
    start_time = time()
    solve_DCGLP(master_env,x̂,t̂, DCGLP_env, sub_env.BSPProblem, sub_env.algo_params.SplitCGLPNormType; time_limit)
    DCGLP_time = time() - start_time

    # add benders cut
    if DCGLP_env.ifsolved == true
        # add split cut
        γ₀, γₓ, γₜ = generate_cut(master_env, DCGLP_env, sub_env.algo_params.StrengthenCutStrategy)
        @constraint(master_env.model, -γ₀ - γₓ'master_env.model[:x] - γₜ*master_env.model[:t][w] >= 0) 
        push!(sub_env.split_info.γ₀s, γ₀)
        push!(sub_env.split_info.γₓs, γₓ)
        push!(sub_env.split_info.γₜs, γₜ)
        generate_cut!(master_env, w,DCGLP_env, sub_env.algo_params.SplitBendersStrategy)
    else
        generate_cut!(master_env, w,DCGLP_env, ALL_SPLIT_BENDERS_STRATEGY)
    end
    
    sub_time,ex = generate_cut(master_env, sub_env, w, ORDINARY_CUTSTRATEGY)
    return DCGLP_time+sub_time,ex
end



function generate_cut!(master_env::AbstractMasterEnv,w,main_env::AbstractDCGLPEnv,::NoSplitBendersStrategy)
end


function generate_cut!(master_env::AbstractMasterEnv,w,main_env::AbstractDCGLPEnv,::AllSplitBendersStrategy)
    for ex1 in main_env.masterconπpoints1
        @constraint(master_env.model,master_env.model[:t][w] >= ex1)
    end
    for ex2 in main_env.masterconπpoints2
        @constraint(master_env.model,master_env.model[:t][w] >= ex2)
    end
end


function generate_cut!(master_env::AbstractMasterEnv,w,main_env::AbstractDCGLPEnv,::TightSplitBendersStrategy)
    λ₁ = dual.(main_env.conπpoints1)
    λ₂ = dual.(main_env.conπpoints2)
    for i in eachindex(λ₁)
        if λ₁[i] > 1e-05
            @constraint(master_env.model, master_env.model[:t][w] >= main_env.masterconπpoints1[i])
        end
    end
    for i in eachindex(λ₂)
        if λ₂[i] > 1e-05
            @constraint(master_env.model, master_env.model[:t][w] >= main_env.masterconπpoints2[i])
        end
    end
end