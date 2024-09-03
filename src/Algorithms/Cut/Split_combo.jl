function generate_cut(
    master_env::AbstractMasterEnv,
    sub_env::AbstractSubEnv,
    ::SplitCutStrategy;
    time_limit = 1000.00)

    # select the index 
    a,b = select_split_set(master_env, sub_env.algo_params.SplitSetSelectionPolicy)

    DCGLP_env = DCGLP(sub_env, a, b, sub_env.algo_params.SplitCGLPNormType)

    # x_out,t_out = master_env.value_x, master_env.value_t
    # x_in,t_in = master_env.solution_in[1], master_env.solution_in[2]

    # x_in_new,t_in_new = (x_out+x_in)/2, (t_out+t_in)/2
    # master_env.solution_in = [x_in_new,t_in_new]
    # x_sep,t_sep = 0.8*x_out+0.2*x_in_new, 0.8*t_out+0.2*t_in_new
    # x̂,t̂ = x_sep,t_sep
    # x_best,t_best = master_env.best_solution[1], master_env.best_solution[2]

    # master_env.best_solution = [0.5*x_best+0.5*_x̂, 0.5*t_best+0.5*_t̂]
    # x̂ = 0.5*x_best+0.5*_x̂
    # t̂ = 0.5*t_best+0.5*_t̂
    # random_number = randn(length(_x̂))
    # # @info random_number
    # _x̂ = _x̂ .+ random_number*0.1
    # _x̂ = min.(_x̂,1)
    # _x̂ = max.(_x̂,0)
    # x̂ = _x̂
    # # t̂ = _t̂
    # t̂ = 0.2*min(sub_env.obj_value,1e10) + 0.8*_t̂
    # @info x̂
    # @info sub_env.obj_value,_t̂,t̂
    # x̂,t̂ = master_env.value_x, master_env.value_t
    start_time = time()

    x̂,t̂ = master_env.value_x, master_env.value_t
    solve_DCGLP(master_env,x̂,t̂, DCGLP_env, sub_env.BSPProblem, sub_env.BSPProblem, sub_env.algo_params.SplitCGLPNormType; time_limit)

    # solve_DCGLP(master_env,x̂,t̂, DCGLP_env, sub_env.BSPProblem, sub_env.BSPProblem2, sub_env.algo_params.SplitCGLPNormType; time_limit)

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


# function generate_cut!(master_env::AbstractMasterEnv,main_env::AbstractDCGLPEnv,::TightSplitBendersStrategy)
#     λ₁ = dual.(main_env.conπpoints1)
#     λ₂ = dual.(main_env.conπpoints2)
#     for i in eachindex(λ₁)
#         if λ₁[i] > 1e-05
#             @constraint(master_env.model, master_env.model[:t] >= main_env.masterconπpoints1[i])
#         end
#     end
#     for i in eachindex(λ₂)
#         if λ₂[i] > 1e-05
#             @constraint(master_env.model, master_env.model[:t] >= main_env.masterconπpoints2[i])
#         end
#     end
# end

#############################################


function select_split_set(master_env::AbstractMasterEnv, ::MostFracIndex)
    gap_x = abs.(master_env.value_x .- 0.5)
    index = findmin(gap_x)[2]
    a = zeros(Int, length(master_env.value_x))
    a[index] = 1
    @info "index = $index"
    return a, 0
end


function select_split_set(master_env::AbstractMasterEnv, ::RandomIndex)
    a = zeros(Int, length(master_env.value_x))
    index = rand(1:length(master_env.value_x))
    a[index] = 1
    @info "index = $index"
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