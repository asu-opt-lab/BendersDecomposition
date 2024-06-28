function solve_DCGLP(
    master_env,
    x̂,
    t̂,
    main_env::AbstractDCGLPEnv, 
    bsp_env::AbstractSubEnv,
    bsp_env2::AbstractSubEnv,
    pConeType::GammaNorm;
    time_limit)

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

    k̂₀s = []
    k̂ₓs = []
    k̂ₜs = []
    v̂₀s = []
    v̂ₓs = []
    v̂ₜs = []
    status1s = []
    status2s = []

    start_time = time()

    # println("#####################solving main problem#####################")
    while true

        k += 1
        main_time_limit = time() - start_time
        set_time_limit_sec(main_env.model, max(time_limit-main_time_limit,1))
        tt = time()
        optimize!(main_env.model)
        @info "main_time = $(time()-tt)"
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

        push!(k̂₀s, k̂₀)
        push!(k̂ₓs, k̂ₓ)
        push!(k̂ₜs, k̂ₜ)
        push!(v̂₀s, v̂₀)
        push!(v̂ₓs, v̂ₓ)
        push!(v̂ₜs, v̂ₜ)

        # @info "k̂ₓ = $k̂ₓ"
        @info "k̂₀ = $k̂₀"
        # @info "v̂ₓ = $v̂ₓ"
        @info "v̂₀ = $v̂₀"
        # @info "k̂ₜ = $k̂ₜ"
        # @info "v̂ₜ = $v̂ₜ"
        # BSP1
        if k̂₀ != 0 || k == 1
        # if k̂₀<=1e-09
                
        
            set_normalized_rhs.(bsp_env.model[:cx], k̂ₓ)
            # set_normalized_rhs.(bsp_env.model[:cx], k̂ₓ./k̂₀)
            set_normalized_rhs.(bsp_env.model[:cb], k̂₀)
        
            # @info "k̂ₓ = $k̂ₓ"
            # @info "k̂₀ = $k̂₀"
            # @info k̂ₓ./k̂₀

            # set_normalized_rhs.(bsp_env.model[:cx], k̂ₓ)
            # set_normalized_rhs.(bsp_env.model[:cb], k̂₀)

            bsp_time_limit = time() - start_time
            # set_time_limit_sec(bsp_env.model, max(time_limit-bsp_time_limit,1))
            tt = time()
            optimize!(bsp_env.model)
            @info "bsp_time1 = $(time()-tt)"
            status1 = dual_status(bsp_env.model)
            # x1 = value.(bsp_env.model[:x])

            if status1 == FEASIBLE_POINT
                g₁ = objective_value(bsp_env.model)
                @info "g₁ = $g₁"
                # ex1 = @expression(main_env.model, g₁ + dual.(bsp_env.model[:cx])⋅(main_env.model[:kₓ]-k̂ₓ) + dual(bsp_env.model[:cb])*(main_env.model[:k₀]-k̂₀) - main_env.model[:kₜ]) 
                ex1 = @expression(main_env.model, dual.(bsp_env.model[:cx])⋅main_env.model[:kₓ] + sum(dual.(bsp_env.model[:cb]))*main_env.model[:k₀] - main_env.model[:kₜ])
                _UB1 = g₁ - k̂ₜ
                # @info "g₁ = $g₁"
                # @info "k̂ₓ = $k̂ₓ"
                # @info "k̂ₜ = $k̂ₜ"
                # @info value.(bsp_env.model[:x])
                # push!(masterconπpoints1, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
                push!(masterconπpoints1, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + sum(dual.(bsp_env.model[:cb]))))

            elseif status1 == INFEASIBILITY_CERTIFICATE 
                @info status1
                g₁ = Inf
                # ex1 = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:kₓ] + dual(bsp_env.model[:cb])*main_env.model[:k₀])
                ex1 = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:kₓ] + sum(dual.(bsp_env.model[:cb]))*main_env.model[:k₀])
                # push!(conπrays1, @expression(master.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
                push!(conπrays1, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + sum(dual.(bsp_env.model[:cb]))))
            else
                g₁ = Inf
                @error "Wrong status1 = $status1"
            end
        else
            # @info "k̂ₓ = $k̂ₓ"
            # @info "k̂₀ = $k̂₀"
            g₁ = 0
        end
        _UB1 = min(_UB1, g₁ - k̂ₜ)

        # BSP2
        if v̂₀ != 0 || k == 1
        # if v̂₀<=1e-09
                # @info "v̂₀ = $v̂₀"
                # @info "v̂ₓ = $v̂ₓ"
            # end
            
            set_normalized_rhs.(bsp_env2.model[:cx], v̂ₓ)
            # set_normalized_rhs.(bsp_env.model[:cx], v̂ₓ./v̂₀)
            set_normalized_rhs.(bsp_env2.model[:cb], v̂₀)
     
      
            # @info "v̂ₓ = $v̂ₓ"
            # @info "v̂₀ = $v̂₀"
            # @info v̂ₓ./v̂₀
            # @info "diff = $(k̂ₓ./k̂₀ - v̂ₓ./v̂₀)"
            bsp_time_limit = time() - start_time
            # set_time_limit_sec(bsp_env2.model, max(time_limit-bsp_time_limit,1))
            tt = time()
            optimize!(bsp_env2.model)
            @info "bsp_time2 = $(time()-tt)"
            status2 = dual_status(bsp_env2.model)
            # x2 = value.(bsp_env.model[:x])
            @info termination_status(bsp_env2.model)
            # @info "diff = $(x1-x2)"
            if status2 == FEASIBLE_POINT && termination_status(bsp_env2.model) !=  TIME_LIMIT
                g₂ = objective_value(bsp_env2.model)
                # @info "g₂ = $g₂"
                # @info "v̂ₓ = $v̂ₓ"
                # @info "v̂ₜ = $v̂ₜ"
                # @info value.(bsp_env.model[:x])
                # ex2 = @expression(main_env.model, g₂ + dual.(bsp_env.model[:cx])⋅(main_env.model[:vₓ]-v̂ₓ) + dual(bsp_env.model[:cb])*(main_env.model[:v₀]-v̂₀) - main_env.model[:vₜ])
                ex2 = @expression(main_env.model, dual.(bsp_env2.model[:cx])⋅main_env.model[:vₓ] + sum(dual.(bsp_env2.model[:cb]))*main_env.model[:v₀] - main_env.model[:vₜ])
                _UB2 = g₂ - v̂ₜ
                # push!(masterconπpoints2, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
                # @constraint(master_env.model, master_env.model[:t] >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
                push!(masterconπpoints2, @expression(master_env.model, dual.(bsp_env2.model[:cx])'master_env.model[:x] + sum(dual.(bsp_env2.model[:cb]))))
            elseif status2 == INFEASIBILITY_CERTIFICATE && termination_status(bsp_env2.model) !=  TIME_LIMIT    
                @info status2      
                g₂ = Inf
                # ex2 = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:vₓ] + dual(bsp_env.model[:cb])*main_env.model[:v₀])
                ex2 = @expression(main_env.model, dual.(bsp_env2.model[:cx])'main_env.model[:vₓ] + sum(dual.(bsp_env2.model[:cb]))*main_env.model[:v₀])
                # push!(conπrays2, @expression(master.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
                # @constraint(master_env.model, 0 >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
                push!(conπrays2, @expression(master_env.model, dual.(bsp_env2.model[:cx])'master_env.model[:x] + sum(dual.(bsp_env2.model[:cb]))))
            else
                g₂ = Inf
                @error "Wrong status2 = $status2"
            end
        else
            # @info "v̂₀ = $v̂₀"
            # @info "v̂ₓ = $v̂ₓ"
            g₂ = 0
        end
        _UB2 = min(_UB2, g₂ - v̂ₜ)
        # push!(status1s, status1)
        # push!(status2s, status2)

        # _UB1 = min(_UB1, g₁ - k̂ₜ)
        # _UB2 = min(_UB2, g₂ - v̂ₜ)
         
        LB = τ̂
        UB = update_UB!(UB,_sx,g₁,g₂,t̂, pConeType)

        @info "Iteration $k: LB = $LB, UB = $UB, _UB1 = $_UB1, _UB2 = $_UB2"

        if ((UB - LB)/abs(UB) <= 1e-3 || (1e-3 >= _UB1 && 1e-3 >= _UB2 )) || (UB - LB) <= 0.01
            main_env.ifsolved = true
            break
        end

        if time_limit <= time() - start_time
            main_env.ifsolved = false
            break
        end

        
        if k̂₀ == 0
            if 1e-3 < _UB1
                push!(conπpoints1, @constraint(main_env.model, 0 >= ex1))
                @info "_add feasible cut 1"
            end
        else
            if status1 == FEASIBLE_POINT 
                if 1e-3 < _UB1
                    push!(conπpoints1, @constraint(main_env.model, 0 >= ex1))
                    @info "add feasible cut 1"
                end
            elseif status1 == INFEASIBILITY_CERTIFICATE
                push!(conπrays1, @constraint(main_env.model, 0 >= 10*ex1))
                @info "add infeasible cut 1"
            end
        end

        if termination_status(bsp_env2.model) !=  TIME_LIMIT
            if v̂₀ == 0
                if 1e-3 < _UB2
                    push!(conπpoints2, @constraint(main_env.model, 0 >= ex2))
                    @info "_add feasible cut 2"
                end
            else
                if status2 == FEASIBLE_POINT
                    if 1e-3 < _UB2
                        push!(conπpoints2, @constraint(main_env.model, 0 >= ex2))
                        @info "add feasible cut 2"
                    end
                elseif status2 == INFEASIBILITY_CERTIFICATE
                    push!(conπrays2, @constraint(main_env.model, 0 >= 10*ex2))
                    @info "add infeasible cut 2"
                end
            end
        end
        

       

    end

    main_env.masterconπpoints1 = masterconπpoints1
    main_env.masterconπpoints2 = masterconπpoints2
    main_env.conπpoints1 = conπpoints1
    main_env.conπpoints2 = conπpoints2

    @info "iteration = $k"
    @info "k̂₀ = $(k̂₀s[end])"
    @info "v̂₀ = $(v̂₀s[end])"
    # @info status1s
    # @info status2s
    point = [k̂ₓs[end]; v̂ₓs[end] ; k̂ₜs[end] ; v̂ₜs[end]; k̂₀s[end]; v̂₀s[end]]
    # for iter in 1:k
    #     # @info "k̂ₓs[$iter] = $(k̂ₓs[iter])"
    #     # @info "v̂ₓs[$iter] = $(v̂ₓs[iter])"
    #     # @info "distance_k_$iter = $(norm(k̂ₓs[iter] - k̂ₓs[end], 2))"
    #     # @info "distance_v_$iter = $(norm(v̂ₓs[iter] - v̂ₓs[end], 2))"
    #     point_iter = [k̂ₓs[iter]; v̂ₓs[iter] ; k̂ₜs[iter] ; v̂ₜs[iter]; k̂₀s[iter]; v̂₀s[iter]]
    #     @info "distance_inf_$iter = $(norm(point_iter - point, Inf))"
    #     # println("distance_L2_$iter = $(norm(point_iter - point, 2))")
    # end
    for iter in 2:k
        # @info "k̂ₓs[$iter] = $(k̂ₓs[iter])"
        # @info "v̂ₓs[$iter] = $(v̂ₓs[iter])"
        # @info "distance_k_$iter = $(norm(k̂ₓs[iter] - k̂ₓs[end], 2))"
        # @info "distance_v_$iter = $(norm(v̂ₓs[iter] - v̂ₓs[end], 2))"
        point_iter_ = [k̂ₓs[iter-1]; v̂ₓs[iter-1] ; k̂ₜs[iter-1] ; v̂ₜs[iter-1]; k̂₀s[iter-1]; v̂₀s[iter-1]]
        point_iter = [k̂ₓs[iter]; v̂ₓs[iter] ; k̂ₜs[iter] ; v̂ₜs[iter]; k̂₀s[iter]; v̂₀s[iter]]
        @info "_distance_L2_$iter = $(norm(point_iter - point_iter_, 2))"
        @info "distance_L2_$iter = $(norm(point_iter - point, 2))"
       # println("distance_inf_$iter = $(norm(point_iter - point, Inf))")
    end
end