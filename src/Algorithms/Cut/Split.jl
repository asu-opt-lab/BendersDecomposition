function generate_cut(
    master_env::AbstractMasterEnv,
    sub_env::CFLPSplitSubEnv,
    ::SplitCutStrategy;
    time_limit = 1000.00)

    # select the index 
    a,b = select_split_set(master_env, sub_env.algo_params.SplitSetSelectionPolicy)

    DCGLP_env = DCGLP(sub_env, a, b, sub_env.algo_params.SplitCGLPNormType)

    start_time = time()
    solve_DCGLP(master_env, DCGLP_env, sub_env.BSPProblem, sub_env.algo_params.SplitCGLPNormType; time_limit)
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
    
    return DCGLP_time+sub_time,ex
end


function solve_DCGLP(
    master_env::AbstractMasterEnv,
    main_env::AbstractDCGLPEnv, 
    bsp_env::AbstractSubEnv,
    ::AbstractNormType;
    time_limit)

    @error "wrong type set"

end


function solve_DCGLP(
    master_env::AbstractMasterEnv,
    main_env::AbstractDCGLPEnv, 
    bsp_env::AbstractSubEnv,
    ::StandardNorm;
    time_limit)

    x̂,t̂ = master_env.value_x, master_env.value_t
    k = 0
    LB = -Inf
    UB = Inf
    conπpoints1 = []
    conπpoints2 = []
    conπrays1 = []
    conπrays2 = []
    _UB1 = Inf
    _UB2 = Inf
    
    set_normalized_rhs.(main_env.model[:conx], x̂)
    set_normalized_rhs.(main_env.model[:cont], t̂)

    start_time = time()
    

    # println("#####################solving main problem#####################")
    while true

        k += 1
        main_time_limit = time() - start_time
        set_time_limit_sec(main_env.model, max(time_limit-main_time_limit,5))
        optimize!(main_env.model)
        k̂₀ = value(main_env.model[:k₀])
        k̂ₓ = value.(main_env.model[:kₓ])
        k̂ₜ = value(main_env.model[:kₜ])
        v̂₀ = value(main_env.model[:v₀])
        v̂ₓ = value.(main_env.model[:vₓ])
        v̂ₜ = value(main_env.model[:vₜ])
        τ̂ = value(main_env.model[:τ])

        # BSP1
        set_normalized_rhs.(bsp_env.model[:cx], k̂ₓ)
        set_normalized_rhs(bsp_env.model[:cb], k̂₀)
        bsp_time_limit = time() - start_time
        set_time_limit_sec(bsp_env.model, max(time_limit-bsp_time_limit,5))
        optimize!(bsp_env.model)
        status1 = dual_status(bsp_env.model)

        if status1 == FEASIBLE_POINT
            g₁ = objective_value(bsp_env.model)
            ex1 = @expression(main_env.model,  -main_env.model[:τ] + g₁ + dual.(bsp_env.model[:cx])⋅(main_env.model[:kₓ]-k̂ₓ) + dual(bsp_env.model[:cb])*(main_env.model[:k₀]-k̂₀) - main_env.model[:kₜ]) 
            _UB1 = g₁ - k̂ₜ
            # @constraint(master_env.model, master_env.model[:t] >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
        elseif status1 == INFEASIBILITY_CERTIFICATE
            g₁ = Inf
            ex1 = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:kₓ] + dual(bsp_env.model[:cb])*main_env.model[:k₀])
            # @constraint(master_env.model, 0 >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
        else
            g₁ = Inf
            @error "Wrong status1: $status1"
        end

        # BSP2
        set_normalized_rhs.(bsp_env.model[:cx], v̂ₓ)
        set_normalized_rhs(bsp_env.model[:cb], v̂₀)
        bsp_time_limit = time() - start_time
        set_time_limit_sec(bsp_env.model, max(time_limit-bsp_time_limit,5))
        optimize!(bsp_env.model)
        status2 = dual_status(bsp_env.model)

        if status2 == FEASIBLE_POINT
            g₂ = objective_value(bsp_env.model)
            ex2 = @expression(main_env.model, -main_env.model[:τ] + g₂ + dual.(bsp_env.model[:cx])⋅(main_env.model[:vₓ]-v̂ₓ) + dual(bsp_env.model[:cb])*(main_env.model[:v₀]-v̂₀) - main_env.model[:vₜ])
            _UB2 = g₂ - v̂ₜ
            # @constraint(master_env.model, master_env.model[:t] >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
        elseif status2 == INFEASIBILITY_CERTIFICATE           
            g₂ = Inf
            ex2 = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:vₓ] + dual(bsp_env.model[:cb])*main_env.model[:v₀])
            # @constraint(master_env.model, 0 >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
        else
            g₂ = Inf
            @error "Wrong status2: $status2"
        end

        _UB1 = min(_UB1, g₁ - k̂ₜ)
        _UB2 = min(_UB2, g₂ - v̂ₜ)


        LB = τ̂
        UB = min(max(_UB1, _UB2),UB)
        @info "Iteration $k: LB = $LB, UB = $UB, _UB1 = $_UB1, _UB2 = $_UB2"

        if (UB - LB) <= 0.01 || (τ̂  >= _UB1 && τ̂  >= _UB2 ) || (UB - LB)/abs(UB) <= 1e-03
            main_env.ifsolved = true
            break
        end

        if time_limit <= time() - start_time
            @info time() - start_time
            main_env.ifsolved = false
            break
        end

        if status1 == FEASIBLE_POINT
            if τ̂  < _UB1
                push!(conπpoints1, @constraint(main_env.model, 0 >= ex1))
            end
        elseif status1 == INFEASIBILITY_CERTIFICATE
            push!(conπrays1, @constraint(main_env.model, 0 >= 10*ex1))
        end

        if status2 == FEASIBLE_POINT
            if τ̂  < _UB2
                push!(conπpoints2, @constraint(main_env.model, 0 >= ex2))
            end
        elseif status2 == INFEASIBILITY_CERTIFICATE
            push!(conπrays2, @constraint(main_env.model, 0 >= 10*ex2))
        end

    end

end


function solve_DCGLP(
    master_env::AbstractMasterEnv,
    main_env::AbstractDCGLPEnv, 
    bsp_env::AbstractSubEnv,
    pConeType::GammaNorm;
    time_limit)

    x̂,t̂ = master_env.value_x, master_env.value_t
    k = 0
    LB = -Inf
    UB = Inf
    conπpoints1 = []
    conπpoints2 = []
    conπrays1 = []
    conπrays2 = []
    masterconπpoints1 = []
    masterconπpoints2 = []
    # masterconπrays1 = []
    # masterconπrays2 = []
    _UB1 = Inf
    _UB2 = Inf
    
    set_normalized_rhs.(main_env.model[:conx], x̂)
    set_normalized_rhs.(main_env.model[:cont], t̂)

    start_time = time()

    # println("#####################solving main problem#####################")
    while true

        k += 1
        main_time_limit = time() - start_time
        set_time_limit_sec(main_env.model, max(time_limit-main_time_limit,1))
        optimize!(main_env.model)
        if termination_status(main_env.model) == TIME_LIMIT
            break
        end
        k̂₀ = value(main_env.model[:k₀])
        k̂ₓ = value.(main_env.model[:kₓ])
        k̂ₜ = value(main_env.model[:kₜ])
        v̂₀ = value(main_env.model[:v₀])
        v̂ₓ = value.(main_env.model[:vₓ])
        v̂ₜ = value(main_env.model[:vₜ])
        τ̂ = value(main_env.model[:τ])
        _sx = value.(main_env.model[:sx])

        # BSP1
        set_normalized_rhs.(bsp_env.model[:cx], k̂ₓ)
        set_normalized_rhs(bsp_env.model[:cb], k̂₀)
        bsp_time_limit = time() - start_time
        set_time_limit_sec(bsp_env.model, max(time_limit-bsp_time_limit,1))
        optimize!(bsp_env.model)
        status1 = dual_status(bsp_env.model)
        
        if status1 == FEASIBLE_POINT
            g₁ = objective_value(bsp_env.model)
            ex1 = @expression(main_env.model, g₁ + dual.(bsp_env.model[:cx])⋅(main_env.model[:kₓ]-k̂ₓ) + dual(bsp_env.model[:cb])*(main_env.model[:k₀]-k̂₀) - main_env.model[:kₜ]) 
            _UB1 = g₁ - k̂ₜ
            push!(masterconπpoints1, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
            # @constraint(master_env.model, master_env.model[:t] >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
        elseif status1 == INFEASIBILITY_CERTIFICATE 
            g₁ = Inf
            ex1 = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:kₓ] + dual(bsp_env.model[:cb])*main_env.model[:k₀])
            # push!(conπrays1, @expression(master.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
            # @constraint(master_env.model, 0 >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
        else
            g₁ = Inf
            @error "Wrong status1"
        end

        # BSP2
        set_normalized_rhs.(bsp_env.model[:cx], v̂ₓ)
        set_normalized_rhs(bsp_env.model[:cb], v̂₀)
        bsp_time_limit = time() - start_time
        set_time_limit_sec(bsp_env.model, max(time_limit-bsp_time_limit,1))
        optimize!(bsp_env.model)
        status2 = dual_status(bsp_env.model)
        
        if status2 == FEASIBLE_POINT
            g₂ = objective_value(bsp_env.model)
            ex2 = @expression(main_env.model, g₂ + dual.(bsp_env.model[:cx])⋅(main_env.model[:vₓ]-v̂ₓ) + dual(bsp_env.model[:cb])*(main_env.model[:v₀]-v̂₀) - main_env.model[:vₜ])
            _UB2 = g₂ - v̂ₜ
            push!(masterconπpoints2, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
            # @constraint(master_env.model, master_env.model[:t] >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
        elseif status2 == INFEASIBILITY_CERTIFICATE           
            g₂ = Inf
            ex2 = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:vₓ] + dual(bsp_env.model[:cb])*main_env.model[:v₀])
            # push!(conπrays2, @expression(master.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
            # @constraint(master_env.model, 0 >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
        else
            g₂ = Inf
            @error "Wrong status2"
        end

        _UB1 = min(_UB1, g₁ - k̂ₜ)
        _UB2 = min(_UB2, g₂ - v̂ₜ)
         
        LB = τ̂
        UB = update_UB!(UB,_sx,g₁,g₂,t̂, pConeType)

        @info "Iteration $k: LB = $LB, UB = $UB, _UB1 = $_UB1, _UB2 = $_UB2"

        if ((UB - LB)/abs(UB) <= 1e-3 || (1e-3 >= _UB1 && 1e-3 >= _UB2 )) 
            main_env.ifsolved = true
            break
        end

        if time_limit <= time() - start_time
            main_env.ifsolved = false
            break
        end
        if status1 == FEASIBLE_POINT
            if 1e-3 < _UB1
                push!(conπpoints1, @constraint(main_env.model, 0 >= ex1))
            end
        elseif status1 == INFEASIBILITY_CERTIFICATE
            push!(conπrays1, @constraint(main_env.model, 0 >= 10*ex1))
        end

        if status2 == FEASIBLE_POINT
            if 1e-3 < _UB2
                push!(conπpoints2, @constraint(main_env.model, 0 >= ex2))
            end
        elseif status2 == INFEASIBILITY_CERTIFICATE
            push!(conπrays2, @constraint(main_env.model, 0 >= 10*ex2))
        end

    end

    main_env.masterconπpoints1 = masterconπpoints1
    main_env.masterconπpoints2 = masterconπpoints2
    main_env.conπpoints1 = conπpoints1
    main_env.conπpoints2 = conπpoints2
end

function update_UB!(UB,_sx,g₁,g₂,t̂,::L1GammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], Inf)) end
function update_UB!(UB,_sx,g₁,g₂,t̂,::L2GammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], 2)) end
function update_UB!(UB,_sx,g₁,g₂,t̂,::LInfGammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], 1)) end


##################################################
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
    end
    for ex2 in main_env.masterconπpoints2
        @constraint(master_env.model,master_env.model[:t] >= ex2)
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