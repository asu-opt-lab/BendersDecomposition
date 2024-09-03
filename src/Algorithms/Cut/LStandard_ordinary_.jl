function solve_DCGLP(
    master_env,
    x̂,
    t̂,
    main_env::AbstractDCGLPEnv, 
    bsp_env::AbstractSubEnv,
    bsp_env2::AbstractSubEnv,
    pConeType::StandardNorm;
    time_limit=100)

    # @info x̂
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

    # k̂₀s = []
    # k̂ₓs = []
    # k̂ₜs = []
    # v̂₀s = []
    # v̂ₓs = []
    # v̂ₜs = []

    start_time = time()

    # println("#####################solving main problem#####################")
    while true

        ##################### DCGLP #####################
        k += 1
        main_time_limit = time() - start_time
        set_time_limit_sec(main_env.model, max(time_limit-main_time_limit,1))
        tt = time()
        optimize!(main_env.model)
        @info "main_time = $(time()-tt)"
        if termination_status(main_env.model) == TIME_LIMIT
            k -= 1
            break
        end
        k̂₀ = value(main_env.model[:k₀])
        k̂ₓ = value.(main_env.model[:kₓ])
        k̂ₜ = value(main_env.model[:kₜ])
        v̂₀ = value(main_env.model[:v₀])
        v̂ₓ = value.(main_env.model[:vₓ])
        v̂ₜ = value(main_env.model[:vₜ])
        τ̂ = value(main_env.model[:τ])

        # push!(k̂₀s, k̂₀)
        # push!(k̂ₓs, k̂ₓ)
        # push!(k̂ₜs, k̂ₜ)
        # push!(v̂₀s, v̂₀)
        # push!(v̂ₓs, v̂ₓ)
        # push!(v̂ₜs, v̂ₜ)

        ##################### BSP1 #####################
        if k̂₀ != 0 #|| k == 1
                
            set_normalized_rhs.(bsp_env.model[:cx], k̂ₓ)
            set_normalized_rhs.(bsp_env.model[:cb], k̂₀)

            bsp_time_limit = time() - start_time
            # set_time_limit_sec(bsp_env.model, max(time_limit-bsp_time_limit,1))
            tt = time()
            optimize!(bsp_env.model)
            @info "bsp_time1 = $(time()-tt)"
            status1 = dual_status(bsp_env.model)
            @info status1
            if status1 == FEASIBLE_POINT && termination_status(bsp_env.model) !=  TIME_LIMIT
                g₁ = objective_value(bsp_env.model)
                ex1 = @expression(main_env.model,  -main_env.model[:τ] + g₁ + dual.(bsp_env.model[:cx])⋅(main_env.model[:kₓ]-k̂ₓ) + dual(bsp_env.model[:cb])*(main_env.model[:k₀]-k̂₀) - main_env.model[:kₜ]) 
                _UB1 = g₁ - k̂ₜ
                push!(masterconπpoints1, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
                @constraint(main_env.model, dual.(bsp_env.model[:cx])⋅main_env.model[:vₓ] + sum(dual.(bsp_env.model[:cb]))*main_env.model[:v₀] - main_env.model[:vₜ] <= 0)

            elseif status1 == INFEASIBILITY_CERTIFICATE && termination_status(bsp_env.model) !=  TIME_LIMIT
                @info status1
                g₁ = Inf
                ex1 = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:kₓ] + dual(bsp_env.model[:cb])*main_env.model[:k₀])
                push!(conπrays1, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + sum(dual.(bsp_env.model[:cb]))))
            else
                g₁ = Inf
                @error "Wrong status1 = $status1"
            end
        else
            g₁ = 0
        end
        _UB1 = min(_UB1, g₁ - k̂ₜ)


        ##################### BSP2 #####################
        if v̂₀ != 0 || k == 1

            set_normalized_rhs.(bsp_env2.model[:cx], v̂ₓ)
            set_normalized_rhs.(bsp_env2.model[:cb], v̂₀)

            bsp_time_limit = time() - start_time
            # set_time_limit_sec(bsp_env2.model, max(time_limit-bsp_time_limit,1))
            tt = time()
            optimize!(bsp_env2.model)
            @info "bsp_time2 = $(time()-tt)"
            status2 = dual_status(bsp_env2.model)
            @info termination_status(bsp_env2.model)

            if status2 == FEASIBLE_POINT && termination_status(bsp_env2.model) !=  TIME_LIMIT
                g₂ = objective_value(bsp_env2.model)
                ex2 = @expression(main_env.model, dual.(bsp_env2.model[:cx])⋅main_env.model[:vₓ] + sum(dual.(bsp_env2.model[:cb]))*main_env.model[:v₀] - main_env.model[:vₜ])
                _UB2 = g₂ - v̂ₜ
                @constraint(main_env.model, dual.(bsp_env2.model[:cx])⋅main_env.model[:kₓ] + sum(dual.(bsp_env2.model[:cb]))*main_env.model[:k₀] - main_env.model[:kₜ] <= 0)
                push!(masterconπpoints2, @expression(master_env.model, dual.(bsp_env2.model[:cx])'master_env.model[:x] + sum(dual.(bsp_env2.model[:cb]))))
            elseif status2 == INFEASIBILITY_CERTIFICATE && termination_status(bsp_env2.model) !=  TIME_LIMIT    
                @info status2      
                g₂ = Inf
                ex2 = @expression(main_env.model, dual.(bsp_env2.model[:cx])'main_env.model[:vₓ] + sum(dual.(bsp_env2.model[:cb]))*main_env.model[:v₀])
                push!(conπrays2, @expression(master_env.model, dual.(bsp_env2.model[:cx])'master_env.model[:x] + sum(dual.(bsp_env2.model[:cb]))))
            else
                g₂ = Inf
                @error "Wrong status2 = $status2"
            end
        else
            g₂ = 0
        end
        _UB2 = min(_UB2, g₂ - v̂ₜ)
         


        ##################### LB and UB #####################
        LB = τ̂
        UB = min(max(_UB1, _UB2),UB)

        @info "Iteration $k: LB = $LB, UB = $UB, _UB1 = $_UB1, _UB2 = $_UB2"


        ##################### check termination #####################
        if ((UB - LB)/abs(UB) <= 1e-3 || (τ̂  >= _UB1 && τ̂  >= _UB2 )) || (UB - LB) <= 0.01 || k >= 30
            main_env.ifsolved = true
            break
        end

        if time_limit <= time() - start_time
            main_env.ifsolved = false
            break
        end

        ##################### add cuts into DCGLP and master problem #####################
        if termination_status(bsp_env.model) !=  TIME_LIMIT
            if k̂₀ == 0
                if τ̂  < _UB1
                    push!(conπpoints1, @constraint(main_env.model, 0 >= ex1))
                    @info "_add feasible cut 1"
                end
            else
                if status1 == FEASIBLE_POINT 
                    if τ̂  < _UB1
                        push!(conπpoints1, @constraint(main_env.model, 0 >= ex1))
                        @info "add feasible cut 1"
                    end
                elseif status1 == INFEASIBILITY_CERTIFICATE
                    push!(conπrays1, @constraint(main_env.model, 0 >= ex1))  # 10*ex1
                    @info "add infeasible cut 1"
                end
            end
        end

        if termination_status(bsp_env2.model) !=  TIME_LIMIT
            if v̂₀ == 0
                if τ̂  < _UB2
                    push!(conπpoints2, @constraint(main_env.model, 0 >= ex2))
                    @info "_add feasible cut 2"
                end
            else
                if status2 == FEASIBLE_POINT
                    if τ̂  < _UB2
                        push!(conπpoints2, @constraint(main_env.model, 0 >= ex2))
                        @info "add feasible cut 2"
                    end
                elseif status2 == INFEASIBILITY_CERTIFICATE
                    push!(conπrays2, @constraint(main_env.model, 0 >= ex2))  # 10*ex2
                    @info "add infeasible cut 2"
                end
            end
        end
        
    end

    ##################### update #####################
    main_env.masterconπpoints1 = masterconπpoints1
    main_env.masterconπpoints2 = masterconπpoints2
    main_env.conπpoints1 = conπpoints1  # calculate dual value later need
    main_env.conπpoints2 = conπpoints2

    @info "iteration = $k"
    @info "k̂₀ = $(k̂₀s[end])"
    @info "v̂₀ = $(v̂₀s[end])"
    optimize!(main_env.model)
    k̂₀ = value(main_env.model[:k₀])
    k̂ₓ = value.(main_env.model[:kₓ])
    k̂ₜ = value(main_env.model[:kₜ])
    v̂₀ = value(main_env.model[:v₀])
    v̂ₓ = value.(main_env.model[:vₓ])
    v̂ₜ = value(main_env.model[:vₜ])
    γₜ = dual(main_env.model[:cont])
    γ₀ = dual(main_env.model[:con0])
    γₓ = dual.(main_env.model[:conx])
    value_κ = -γₜ*k̂ₜ/k̂₀ - γₓ'k̂ₓ/k̂₀ - γ₀
    value_ν = -γₜ*v̂ₜ/v̂₀ - γₓ'v̂ₓ/v̂₀ - γ₀
    @info "value_κ = $value_κ"
    @info "value_ν = $value_ν"
end








# @info status1s
    # @info status2s
    # point = [k̂ₓs[end]; v̂ₓs[end] ; k̂ₜs[end] ; v̂ₜs[end]; k̂₀s[end]; v̂₀s[end]]
    # # for iter in 1:k
    # #     # @info "k̂ₓs[$iter] = $(k̂ₓs[iter])"
    # #     # @info "v̂ₓs[$iter] = $(v̂ₓs[iter])"
    # #     # @info "distance_k_$iter = $(norm(k̂ₓs[iter] - k̂ₓs[end], 2))"
    # #     # @info "distance_v_$iter = $(norm(v̂ₓs[iter] - v̂ₓs[end], 2))"
    # #     point_iter = [k̂ₓs[iter]; v̂ₓs[iter] ; k̂ₜs[iter] ; v̂ₜs[iter]; k̂₀s[iter]; v̂₀s[iter]]
    # #     @info "distance_inf_$iter = $(norm(point_iter - point, Inf))"
    # #     # println("distance_L2_$iter = $(norm(point_iter - point, 2))")
    # # end
    # for iter in 2:k
    #     # @info "k̂ₓs[$iter] = $(k̂ₓs[iter])"
    #     # @info "v̂ₓs[$iter] = $(v̂ₓs[iter])"
    #     # @info "distance_k_$iter = $(norm(k̂ₓs[iter] - k̂ₓs[end], 2))"
    #     # @info "distance_v_$iter = $(norm(v̂ₓs[iter] - v̂ₓs[end], 2))"
    #     point_iter_ = [k̂ₓs[iter-1]; v̂ₓs[iter-1] ; k̂ₜs[iter-1] ; v̂ₜs[iter-1]; k̂₀s[iter-1]; v̂₀s[iter-1]]
    #     point_iter = [k̂ₓs[iter]; v̂ₓs[iter] ; k̂ₜs[iter] ; v̂ₜs[iter]; k̂₀s[iter]; v̂₀s[iter]]
    #     # @info "_distance_L2_$iter = $(norm(point_iter - point_iter_, 2))"
    #     # @info "distance_L2_$iter = $(norm(point_iter - point, 2))"
    #    # println("distance_inf_$iter = $(norm(point_iter - point, Inf))")
    # end