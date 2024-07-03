function solve_DCGLP(
    master_env,
    x̂,
    t̂,
    main_env::AbstractDCGLPEnv, 
    bsp_env::AbstractSubEnv,
    ::StandardNorm;
    time_limit)

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
    k̂₀s = []
    k̂ₓs = []
    k̂ₜs = []
    v̂₀s = []
    v̂ₓs = []
    v̂ₜs = []


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

        push!(k̂₀s, k̂₀)
        push!(k̂ₓs, k̂ₓ)
        push!(k̂ₜs, k̂ₜ)
        push!(v̂₀s, v̂₀)
        push!(v̂ₓs, v̂ₓ)
        push!(v̂ₜs, v̂ₜ)


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
        # @info "Iteration $k: LB = $LB, UB = $UB, _UB1 = $_UB1, _UB2 = $_UB2"

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