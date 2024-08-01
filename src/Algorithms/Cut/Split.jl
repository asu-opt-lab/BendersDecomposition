function generate_cut(
    master_env::AbstractMasterEnv,
    sub_env::AbstractSubEnv,
    ::SplitCutStrategy;
    time_limit = 1000.00)

    # select the index 
    a,b = select_split_set(master_env, sub_env.algo_params.SplitSetSelectionPolicy)

    DCGLP_env = DCGLP(sub_env, a, b, sub_env.algo_params.SplitCGLPNormType)

    x̂,t̂ = master_env.value_x, master_env.value_t
    start_time = time()
    # BSPProblem = generate_BSPProblem(sub_env.data, solver=:Gurobi)
    # BSPProblem2 = generate_BSPProblem(sub_env.data, solver=:Gurobi)
    # BSPProblem = generate_BSPProblem_Advanced(sub_env.data)
    # BSPProblem2 = generate_BSPProblem_Advanced(sub_env.data)

    # solve_DCGLP(master_env,x̂,t̂, DCGLP_env, sub_env.BSPProblem, sub_env.algo_params.SplitCGLPNormType; time_limit)
    # solve_DCGLP(master_env,x̂,t̂, DCGLP_env, sub_env.BSPProblem, sub_env.BSPProblem2, sub_env.algo_params.SplitCGLPNormType; time_limit)
    # solve_DCGLP(master_env,x̂,t̂, DCGLP_env, BSPProblem, BSPProblem2, sub_env.algo_params.SplitCGLPNormType; time_limit)

    solve_DCGLP(master_env,x̂,t̂, DCGLP_env, sub_env.BSPProblem, sub_env.BSPProblem2, sub_env.algo_params.SplitCGLPNormType; time_limit)
    # solve_DCGLP(master_env,x̂,t̂, DCGLP_env, sub_env, sub_env.algo_params.SplitCGLPNormType; time_limit)

    DCGLP_time = time() - start_time

    # add benders cut
    if DCGLP_env.ifsolved == true
        # add split cut
        γ₀, γₓ, γₜ = generate_cut(master_env, DCGLP_env, sub_env.algo_params.StrengthenCutStrategy)
        @constraint(master_env.model, -γ₀ - γₓ'master_env.model[:x] - γₜ*master_env.model[:t] >= 0) 
        push!(sub_env.split_info.γ₀s, γ₀)
        push!(sub_env.split_info.γₓs, γₓ)
        push!(sub_env.split_info.γₜs, γₜ)
        generate_cut!(master_env, DCGLP_env, sub_env.algo_params.SplitBendersStrategy)
    else
        generate_cut!(master_env, DCGLP_env, ALL_SPLIT_BENDERS_STRATEGY)
    end
    
    sub_time,ex = generate_cut(master_env, sub_env, ORDINARY_CUTSTRATEGY)
    # sub_time,ex = 0,0
    # sub_env.obj_value = Inf
    return DCGLP_time+sub_time,ex
end


# function solve_DCGLP(
#     master_env,
#     x̂,
#     t̂,
#     main_env::AbstractDCGLPEnv, 
#     bsp_env::AbstractSubEnv,
#     ::AbstractNormType;
#     time_limit)

#     @error "wrong type set"

# end







####################################################
#### strengthening procedure
function generate_cut(master_env::AbstractMasterEnv,main_env::AbstractDCGLPEnv,::AbstractSplitBendersPolicy)
    @error "wrong type set"
end


function generate_cut(master_env::AbstractMasterEnv,main_env::AbstractDCGLPEnv,::SplitPureCutStrategy)
    optimize!(main_env.model)
    γₜ = dual(main_env.model[:cont])
    γ₀ = dual(main_env.model[:con0])
    γₓ = dual.(main_env.model[:conx])
    # @constraint(master_env.model, -γ₀ - γₓ'master_env.model[:x] - γₜ*master_env.model[:t] >= 0) 
    return γ₀, γₓ, γₜ
end

function generate_cut(master_env::AbstractMasterEnv,main_env::AbstractDCGLPEnv,::SplitStrengthenCutStrategy)
    optimize!(main_env.model)
    σ₁ = dual.(main_env.model[:consigma1])
    σ₂ = dual.(main_env.model[:consigma2])
    if σ₁ == 0 && σ₂ == 0
        γₜ = dual(main_env.model[:cont])
        γ₀ = dual(main_env.model[:con0])
        γₓ = dual.(main_env.model[:conx])
    else
        γₓ = dual.(main_env.model[:conx]) 
        γ₁ = γₓ .- dual.(main_env.model[:conv1])
        γ₂ = γₓ .- dual.(main_env.model[:conv2])
        m = (γ₁ .- γ₂) / (σ₂ + σ₁)
        m_lb = floor.(m)
        m_ub = ceil.(m)
        γₓ = min.(γ₁-σ₁*m_lb, γ₂+σ₂*m_ub)
        γₜ = dual(main_env.model[:cont])
        γ₀ = dual(main_env.model[:con0])
    end 
    return γ₀, γₓ, γₜ
    # @constraint(master_env.model, -γ₀ - γₓ'master_env.model[:x] - γₜ*master_env.model[:t] >= 0) 
end


#############################################

function generate_cut!(master_env::AbstractMasterEnv,main_env::AbstractDCGLPEnv,::AbstractSplitBendersPolicy)
    @error "wrong type set"
end

function generate_cut!(master_env::AbstractMasterEnv,main_env::AbstractDCGLPEnv,::NoSplitBendersStrategy)
end


function generate_cut!(master_env::AbstractMasterEnv,main_env::AbstractDCGLPEnv,::AllSplitBendersStrategy)
    for ex1 in main_env.masterconπpoints1
        @constraint(master_env.model,master_env.model[:t] >= ex1)
        # @constraint(master_env.model,0 >= ex1)
        # @info ex1
    end
    for ex2 in main_env.masterconπpoints2
        @constraint(master_env.model,master_env.model[:t] >= ex2)
        # @constraint(master_env.model,0 >= ex2)
        # @info ex2
    end
end


function generate_cut!(master_env::AbstractMasterEnv,main_env::AbstractDCGLPEnv,::TightSplitBendersStrategy)
    λ₁ = dual.(main_env.conπpoints1)
    λ₂ = dual.(main_env.conπpoints2)
    for i in eachindex(λ₁)
        if λ₁[i] > 1e-05
            @constraint(master_env.model, master_env.model[:t] >= main_env.masterconπpoints1[i])
        end
    end
    for i in eachindex(λ₂)
        if λ₂[i] > 1e-05
            @constraint(master_env.model, master_env.model[:t] >= main_env.masterconπpoints2[i])
        end
    end
end

#############################################


function select_split_set(master_env::AbstractMasterEnv, ::MostFracIndex)
    gap_x = abs.(master_env.value_x .- 0.5)
    index = findmin(gap_x)[2]
    a = zeros(Int, length(master_env.value_x))
    a[index] = 1
    return a, 0
end


function select_split_set(master_env::AbstractMasterEnv, ::RandomIndex)
    a = zeros(Int, length(master_env.value_x))
    a[rand(1:length(master_env.value_x))] = 1
    return a, 0
end


# function select_split_set(master_env::AbstractMasterEnv, split_info, ::MostFracIndex)
#     gap_x = abs.(master_env.value_x .- 0.5)
#     index = findmin(gap_x)[2]
#     @info index split_info.indices
#     @info split_info.ifaddall
#     if split_info.indices == []
#         push!(split_info.indices, index)
#         return index
#     end
#     if index >= findmax(split_info.indices)[1]
#         push!(split_info.indices, index)
#         split_info.ifaddall = true
#     end
#     return index
# end