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
        l = length(k̂ₓ)
            # set_normalized_rhs.(bsp_env.model[:cx], k̂ₓ)
            for i in eachindex(k̂ₓ)
                set_normalized_rhs(bsp_env.cconstr[i], k̂ₓ[i])
                set_normalized_rhs(bsp_env.cconstr[l+i], -k̂ₓ[i])
            end
            # set_normalized_rhs.(bsp_env.model[:cx], k̂ₓ./k̂₀)
            set_normalized_rhs.(bsp_env.model[:cb], k̂₀)
            set_normalized_rhs(bsp_env.oconstr, -k̂ₜ)
            # @info "k̂ₓ = $k̂ₓ"
            # @info "k̂₀ = $k̂₀"
            # @info k̂ₓ./k̂₀

            # set_normalized_rhs.(bsp_env.model[:cx], k̂ₓ)
            # set_normalized_rhs.(bsp_env.model[:cb], k̂₀)

            bsp_time_limit = time() - start_time
            set_time_limit_sec(bsp_env.model, max(time_limit-bsp_time_limit,1))
            tt = time()
            optimize!(bsp_env.model)
            @info "bsp_time1 = $(time()-tt)"
            status1 = dual_status(bsp_env.model)
            # x1 = value.(bsp_env.model[:x])

            if status1 == FEASIBLE_POINT
                subObjVal = JuMP.objective_value(bsp_env.model) 
                ex1 = @expression(main_env.model,  
                subObjVal 
                - dual(bsp_env.oconstr)*(main_env.model[:kₜ] - k̂ₜ)
                - sum(dual(bsp_env.cconstr[i]) * (main_env.model[:kₓ][i] - k̂ₓ[i]) for i in eachindex(k̂ₓ)) 
                - sum(dual(bsp_env.cconstr[l+i]) * (main_env.model[:kₓ][i] - k̂ₓ[i]) for i in eachindex(k̂ₓ))
                + sum(dual.(bsp_env.model[:cb]))*(main_env.model[:k₀]-k̂₀))   
            else
                @error "dual of sub is neither feasible nor infeasible certificate: $status"
                throw(-1)
            end
       
        if abs(subObjVal) <= 1e-06
            g₁ = value(bsp_env.obj)#sum(data.costs[i,j] * data.demands[j] * JuMP.value(sub_env.model[:y][i,j]) for i in 1:data.n_facilities, j in 1:data.n_customers)  
        else 
            g₁ = Inf
        end
        _UB1 = min(_UB1, g₁ - k̂ₜ)

        # BSP2
        l = length(v̂ₓ)
            # set_normalized_rhs.(bsp_env2.model[:cx], v̂ₓ)
            for i in eachindex(v̂ₓ)
                set_normalized_rhs(bsp_env2.cconstr[i], v̂ₓ[i])
                set_normalized_rhs(bsp_env2.cconstr[l+i], -v̂ₓ[i])
            end
            # set_normalized_rhs.(bsp_env2.model[:cx], v̂ₓ./v̂₀)
            set_normalized_rhs.(bsp_env2.model[:cb], v̂₀)
            set_normalized_rhs(bsp_env2.oconstr, -v̂ₜ)
            # @info "v̂ₓ = $v̂ₓ"
            # @info "v̂₀ = $v̂₀"
            # @info v̂ₓ./v̂₀

            # set_normalized_rhs.(bsp_env2.model[:cx], v̂ₓ)
            # set_normalized_rhs.(bsp_env2.model[:cb], v̂₀)

            bsp_time_limit = time() - start_time
            set_time_limit_sec(bsp_env2.model, max(time_limit-bsp_time_limit,1))
            tt = time()
            optimize!(bsp_env2.model)
            @info "bsp_time2 = $(time()-tt)"
            status2 = dual_status(bsp_env2.model)
            # x2 = value.(bsp_env2.model[:x])

            if status2 == FEASIBLE_POINT
                subObjVal = JuMP.objective_value(bsp_env2.model) 
                ex2 = @expression(main_env.model,  
                subObjVal 
                - dual(bsp_env2.oconstr)*(main_env.model[:vₜ] - v̂ₜ)
                - sum(dual(bsp_env2.cconstr[i]) * (main_env.model[:vₓ][i] - v̂ₓ[i]) for i in eachindex(v̂ₓ)) 
                - sum(dual(bsp_env2.cconstr[l+i]) * (main_env.model[:vₓ][i] - v̂ₓ[i]) for i in eachindex(v̂ₓ))
                + sum(dual.(bsp_env2.model[:cb]))*(main_env.model[:v₀]-v̂₀))   
            else
                @error "dual of sub is neither feasible nor infeasible certificate: $status"
                throw(-1)
            end

        if abs(subObjVal) <= 1e-06
            g₂ = value(bsp_env2.obj)#sum(data.costs[i,j] * data.demands[j] * JuMP.value(sub_env.model[:y][i,j]) for i in 1:data.n_facilities, j in 1:data.n_customers)  
        else 
            g₂ = Inf
        end

        _UB2 = min(_UB2, g₂ - v̂ₜ)
         
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

        
        # if k̂₀ == 0
        #     if 1e-3 < _UB1
        #         push!(conπpoints1, @constraint(main_env.model, 0 >= ex1))
        #         @info "_add feasible cut 1"
        #     end
        # else
            if status1 == FEASIBLE_POINT 
                if 1e-3 < _UB1
                    push!(conπpoints1, @constraint(main_env.model, 0 >= ex1))
                    @info "add feasible cut 1"
                end
            elseif status1 == INFEASIBILITY_CERTIFICATE
                push!(conπrays1, @constraint(main_env.model, 0 >= 10*ex1))
                @info "add infeasible cut 1"
            end
        # end

        # if v̂₀ == 0
        #     if 1e-3 < _UB2
        #         push!(conπpoints2, @constraint(main_env.model, 0 >= ex2))
        #         @info "_add feasible cut 2"
        #     end
        # else
            if status2 == FEASIBLE_POINT
                if 1e-3 < _UB2
                    push!(conπpoints2, @constraint(main_env.model, 0 >= ex2))
                    @info "add feasible cut 2"
                end
            elseif status2 == INFEASIBILITY_CERTIFICATE
                push!(conπrays2, @constraint(main_env.model, 0 >= 10*ex2))
                @info "add infeasible cut 2"
            end
        # end
        

       

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
    # point = [k̂ₓs[end]; v̂ₓs[end] ; k̂ₜs[end] ; v̂ₜs[end]; k̂₀s[end]; v̂₀s[end]]
    # for iter in 1:k
    #     # @info "k̂ₓs[$iter] = $(k̂ₓs[iter])"
    #     # @info "v̂ₓs[$iter] = $(v̂ₓs[iter])"
    #     # @info "distance_k_$iter = $(norm(k̂ₓs[iter] - k̂ₓs[end], 2))"
    #     # @info "distance_v_$iter = $(norm(v̂ₓs[iter] - v̂ₓs[end], 2))"
    #     point_iter = [k̂ₓs[iter]; v̂ₓs[iter] ; k̂ₜs[iter] ; v̂ₜs[iter]; k̂₀s[iter]; v̂₀s[iter]]
    #     @info "distance_inf_$iter = $(norm(point_iter - point, Inf))"
    #     # println("distance_L2_$iter = $(norm(point_iter - point, 2))")
    # end
    # for iter in 1:k
    #     # @info "k̂ₓs[$iter] = $(k̂ₓs[iter])"
    #     # @info "v̂ₓs[$iter] = $(v̂ₓs[iter])"
    #     # @info "distance_k_$iter = $(norm(k̂ₓs[iter] - k̂ₓs[end], 2))"
    #     # @info "distance_v_$iter = $(norm(v̂ₓs[iter] - v̂ₓs[end], 2))"
    #     point_iter = [k̂ₓs[iter]; v̂ₓs[iter] ; k̂ₜs[iter] ; v̂ₜs[iter]; k̂₀s[iter]; v̂₀s[iter]]
    #     @info "distance_L2_$iter = $(norm(point_iter - point, 2))"
    #    # println("distance_inf_$iter = $(norm(point_iter - point, Inf))")
    # end
end