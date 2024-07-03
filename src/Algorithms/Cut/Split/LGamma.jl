function solve_DCGLP(
    master_env,
    x̂,
    t̂,
    main_env::AbstractDCGLPEnv, 
    bsp_env::AbstractSubEnv,
    pConeType::GammaNorm;
    time_limit)

    # @info "x̂ = $x̂"

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

        push!(k̂₀s, k̂₀)
        push!(k̂ₓs, k̂ₓ)
        push!(k̂ₜs, k̂ₜ)
        push!(v̂₀s, v̂₀)
        push!(v̂ₓs, v̂ₓ)
        push!(v̂ₜs, v̂ₜ)

        @info "k̂₀ = $k̂₀"
        @info "v̂₀ = $v̂₀"

        # BSP1
        set_normalized_rhs.(bsp_env.model[:cx], k̂ₓ)
        set_normalized_rhs.(bsp_env.model[:cb], k̂₀)


        # set_normalized_rhs.(bsp_env.model[:cx], k̂ₓ)
        # set_normalized_rhs.(bsp_env.model[:cb], k̂₀)

        bsp_time_limit = time() - start_time
        set_time_limit_sec(bsp_env.model, max(time_limit-bsp_time_limit,1))
        optimize!(bsp_env.model)
        status1 = dual_status(bsp_env.model)
        
        if status1 == FEASIBLE_POINT
            g₁ = objective_value(bsp_env.model)
            # ex1 = @expression(main_env.model, g₁ + dual.(bsp_env.model[:cx])⋅(main_env.model[:kₓ]-k̂ₓ) + dual(bsp_env.model[:cb])*(main_env.model[:k₀]-k̂₀) - main_env.model[:kₜ]) 
            ex1 = @expression(main_env.model, dual.(bsp_env.model[:cx])⋅main_env.model[:kₓ] + sum(dual.(bsp_env.model[:cb]))*main_env.model[:k₀] - main_env.model[:kₜ])
            _UB1 = g₁ - k̂ₜ
            # push!(masterconπpoints1, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
            push!(masterconπpoints1, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + sum(dual.(bsp_env.model[:cb]))))

        elseif status1 == INFEASIBILITY_CERTIFICATE 
            g₁ = Inf
            # ex1 = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:kₓ] + dual(bsp_env.model[:cb])*main_env.model[:k₀])
            ex = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:kₓ] + sum(dual.(bsp_env.model[:cb]))*main_env.model[:k₀])
            # push!(conπrays1, @expression(master.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
            push!(conπrays1, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + sum(dual.(bsp_env.model[:cb]))))
        else
            g₁ = Inf
            @error "Wrong status1 = $status1"
        end

        # BSP2
        set_optimizer_attribute(bsp_env.model, "Method", 2)
        # @info "v̂ₓ = $v̂ₓ"
        # @info "v̂₀ = $v̂₀"
        new_value = v̂ₓ ./ v̂₀
        # @info "diff = $(new_value - x̂)"
        set_normalized_rhs.(bsp_env.model[:cx], v̂ₓ)
        set_normalized_rhs.(bsp_env.model[:cb], v̂₀)
        bsp_time_limit = time() - start_time
        set_time_limit_sec(bsp_env.model, max(time_limit-bsp_time_limit,1))
        optimize!(bsp_env.model)
        status2 = dual_status(bsp_env.model)
        
        if status2 == FEASIBLE_POINT
            g₂ = objective_value(bsp_env.model)
            # ex2 = @expression(main_env.model, g₂ + dual.(bsp_env.model[:cx])⋅(main_env.model[:vₓ]-v̂ₓ) + dual(bsp_env.model[:cb])*(main_env.model[:v₀]-v̂₀) - main_env.model[:vₜ])
            ex2 = @expression(main_env.model, dual.(bsp_env.model[:cx])⋅main_env.model[:vₓ] + sum(dual.(bsp_env.model[:cb]))*main_env.model[:v₀] - main_env.model[:vₜ])
            _UB2 = g₂ - v̂ₜ
            # push!(masterconπpoints2, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
            # @constraint(master_env.model, master_env.model[:t] >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
            push!(masterconπpoints2, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + sum(dual.(bsp_env.model[:cb]))))
        elseif status2 == INFEASIBILITY_CERTIFICATE           
            g₂ = Inf
            # ex2 = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:vₓ] + dual(bsp_env.model[:cb])*main_env.model[:v₀])
            ex2 = @expression(main_env.model, dual.(bsp_env.model[:cx])'main_env.model[:vₓ] + sum(dual.(bsp_env.model[:cb]))*main_env.model[:v₀])
            # push!(conπrays2, @expression(master.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb])))
            # @constraint(master_env.model, 0 >= dual.(bsp_env.model[:cx])'master_env.model[:x] + dual(bsp_env.model[:cb]))
            push!(conπrays2, @expression(master_env.model, dual.(bsp_env.model[:cx])'master_env.model[:x] + sum(dual.(bsp_env.model[:cb]))))
        else
            g₂ = Inf
            @error "Wrong status2 = $status2"
        end

        push!(status1s, status1)
        push!(status2s, status2)

        _UB1 = min(_UB1, g₁ - k̂ₜ)
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

    # @info status1s
    # @info status2s
    point = [k̂ₓs[end]; v̂ₓs[end] ; k̂ₜs[end] ; v̂ₜs[end]; k̂₀s[end]; v̂₀s[end]]
    for iter in 1:k
        # @info "k̂ₓs[$iter] = $(k̂ₓs[iter])"
        # @info "v̂ₓs[$iter] = $(v̂ₓs[iter])"
        # @info "distance_k_$iter = $(norm(k̂ₓs[iter] - k̂ₓs[end], 2))"
        # @info "distance_v_$iter = $(norm(v̂ₓs[iter] - v̂ₓs[end], 2))"
        point_iter = [k̂ₓs[iter]; v̂ₓs[iter] ; k̂ₜs[iter] ; v̂ₜs[iter]; k̂₀s[iter]; v̂₀s[iter]]
        @info "distance_inf_$iter = $(norm(point_iter - point, Inf))"
        # println("distance_L2_$iter = $(norm(point_iter - point, 2))")
    end
    for iter in 1:k
        # @info "k̂ₓs[$iter] = $(k̂ₓs[iter])"
        # @info "v̂ₓs[$iter] = $(v̂ₓs[iter])"
        # @info "distance_k_$iter = $(norm(k̂ₓs[iter] - k̂ₓs[end], 2))"
        # @info "distance_v_$iter = $(norm(v̂ₓs[iter] - v̂ₓs[end], 2))"
        point_iter = [k̂ₓs[iter]; v̂ₓs[iter] ; k̂ₜs[iter] ; v̂ₜs[iter]; k̂₀s[iter]; v̂₀s[iter]]
        @info "distance_L2_$iter = $(norm(point_iter - point, 2))"
        # println("distance_inf_$iter = $(norm(point_iter - point, Inf))")
    end
end

function update_UB!(UB,_sx,g₁,g₂,t̂,::L1GammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], Inf)) end
function update_UB!(UB,_sx,g₁,g₂,t̂,::L2GammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], 2)) end
function update_UB!(UB,_sx,g₁,g₂,t̂,::LInfGammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], 1)) end
