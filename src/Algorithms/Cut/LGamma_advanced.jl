function solve_DCGLP(
    master_env,
    x̂,
    t̂,
    main_env::AbstractDCGLPEnv, 
    bsp_env1::AbstractSubEnv,
    bsp_env2::AbstractSubEnv,
    pConeType::GammaNorm;
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
    masterconπrays1 = []
    masterconπrays2 = []
    _UB1 = Inf
    _UB2 = Inf
    
    
    set_normalized_rhs.(main_env.model[:conx], x̂)
    set_normalized_rhs.(main_env.model[:cont], t̂)
    coeff = master_env.obj_value
    # @constraint(main_env.model, conk, master_env.data.fixed_costs'main_env.model[:kₓ] + main_env.model[:kₜ] >= coeff*main_env.model[:k₀])
    # @constraint(main_env.model, conv, master_env.data.fixed_costs'main_env.model[:vₓ] + main_env.model[:vₜ] >= coeff*main_env.model[:v₀])
    # @constraint(main_env.model, conk, main_env.model[:kₜ] >= t̂*main_env.model[:k₀])
    # @constraint(main_env.model, conv, main_env.model[:vₜ] >= t̂*main_env.model[:v₀])
   
    # k̂₀s = []
    # k̂ₓs = []
    # k̂ₜs = []
    # v̂₀s = []
    # v̂ₓs = []
    # v̂ₜs = []

    start_time = time()

    # println("#####################solving main problem#####################")
    while true
        # if k>=1
        #     if 1<=k<=50
        #         delta = 1000
        #     elseif 51<=k<=100
        #         delta = 100
        #     elseif 101<=k<=150
        #         delta = 10
        #     end
        #     @constraint(main_env.model, trust_region, [delta; main_env.model[:kₓ] .- k̂ₓ; main_env.model[:vₓ].-v̂ₓ; main_env.model[:k₀]- k̂₀;main_env.model[:v₀]-v̂₀] in MOI.NormInfinityCone(3 + 2*length(x̂)))
        # end
        ##################### DCGLP #####################
        k += 1
        main_time_limit = time() - start_time
        # set_time_limit_sec(main_env.model, max(time_limit-main_time_limit,1))
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
        _sx = value.(main_env.model[:sx])
        # @info main_env.model[:trust_region]
        
        @info "k̂₀ = $k̂₀"
        @info "v̂₀ = $v̂₀"

        # push!(k̂₀s, k̂₀)
        # push!(k̂ₓs, k̂ₓ)
        # push!(k̂ₜs, k̂ₜ)
        # push!(v̂₀s, v̂₀)
        # push!(v̂ₓs, v̂ₓ)
        # push!(v̂ₜs, v̂ₜ)

        ##################### BSP1 #####################
        if k̂₀ != 0 #|| k == 1
            # @info "k̂ₜ/k̂₀ = $(k̂ₜ/k̂₀)"
            set_normalized_rhs.(bsp_env1.model[:cb], k̂₀)
            set_normalized_rhs.(bsp_env1.model[:cbb], -k̂₀)
            set_normalized_rhs.(bsp_env1.cconstr, k̂ₓ)
            set_normalized_rhs(bsp_env1.oconstr, -k̂ₜ)

            bsp_time_limit = time() - start_time
            # set_time_limit_sec(bsp_env.model, max(time_limit-bsp_time_limit,1))
            tt = time()
            optimize!(bsp_env1.model)
            @info "bsp_time1 = $(time()-tt)"
            
            @info status1 = dual_status(bsp_env1.model)
            # if status1 == FEASIBLE_POINT
                g₁ = objective_value(bsp_env1.model)
                @info dual(bsp_env1.oconstr)
                push!(conπpoints1, [dual.(bsp_env1.model[:cb]), dual.(bsp_env1.model[:cbb]), dual.(bsp_env1.cconstr), dual(bsp_env1.oconstr)])
                # ex1 = @expression(main_env.model, 
                # dual.(bsp_env1.sub_constr)'bsp_env1.sub_rhs*main_env.model[:k₀]
                # + dual.(bsp_env1.cconstr)'main_env.model[:kₓ]
                # - main_env.model[:kₜ])
                # push!(masterconπpoints1, @expression(master_env.model, dual.(bsp_env1.model[:bcon])+dual.(bsp_env1.cconstr)'master_env.model[:x] - dual(bsp_env1.oconstr)*master_env.model[:t]))

            # elseif status1 == INFEASIBILITY_CERTIFICATE 
            #     @info status1
            #     g₁ = Inf
            #     push!(conπrays1, [dual.(bsp_env1.sub_constr)'bsp_env1.sub_rhs,dual.(bsp_env1.cconstr)])
            #     # ex1 = @expression(main_env.model, 
            #     # dual.(bsp_env1.sub_constr)'bsp_env1.sub_rhs*main_env.model[:k₀]
            #     # + dual.(bsp_env1.cconstr)'main_env.model[:kₓ])
            #     # push!(conπrays1, @expression(master_env.model, dual.(bsp_env1.sub_constr)'bsp_env1.sub_rhs + dual.(bsp_env1.cconstr)'master_env.model[:x]))
            # else
            #     g₁ = Inf
            #     @error "Wrong status1 = $status1"
            # end
        else
            g₁ = 0
        end
        # _UB1 = min(_UB1, k̂₀*g₁ - k̂ₜ)
        _UB1 = g₁ 

        ##################### BSP2 #####################
        if v̂₀ != 0 #|| k == 1
            # @info "v̂ₜ/v̂₀ = $(v̂ₜ/v̂₀)"
            set_normalized_rhs.(bsp_env2.model[:cb], v̂₀)
            set_normalized_rhs.(bsp_env2.model[:cbb], -v̂₀)
            set_normalized_rhs.(bsp_env2.cconstr, v̂ₓ)
            set_normalized_rhs(bsp_env2.oconstr, -v̂ₜ)

            bsp_time_limit = time() - start_time
            # set_time_limit_sec(bsp_env.model, max(time_limit-bsp_time_limit,1))
            tt = time()
            optimize!(bsp_env2.model)
            @info "bsp_time2 = $(time()-tt)"
            @info status2 = dual_status(bsp_env2.model)

            # if status2 == FEASIBLE_POINT
                g₂ = objective_value(bsp_env2.model)
                push!(conπpoints2, [dual.(bsp_env2.model[:cb]), dual.(bsp_env2.model[:cbb]), dual.(bsp_env2.cconstr), dual(bsp_env2.oconstr)])
                # ex2 = @expression(main_env.model, 
                # dual.(bsp_env2.sub_constr)'bsp_env2.sub_rhs*main_env.model[:v₀]
                # + dual.(bsp_env2.cconstr)'main_env.model[:vₓ]
                # - main_env.model[:vₜ])
                # push!(masterconπpoints2, @expression(master_env.model, dual.(bsp_env2.model[:bcon])+dual.(bsp_env2.cconstr)'master_env.model[:x] - dual(bsp_env2.oconstr)*master_env.model[:t]))
        
            # elseif status2 == INFEASIBILITY_CERTIFICATE && termination_status(bsp_env2.model) !=  TIME_LIMIT
            #     g₂ = Inf
            #     push!(conπrays2, [dual.(bsp_env2.sub_constr)'bsp_env2.sub_rhs,dual.(bsp_env2.cconstr)])
            #     # ex2 = @expression(main_env.model, 
            #     # dual.(bsp_env2.sub_constr)'bsp_env2.sub_rhs*main_env.model[:v₀]
            #     # + dual.(bsp_env2.cconstr)'main_env.model[:vₓ])
            #     # push!(conπrays2, @expression(master_env.model, dual.(bsp_env2.sub_constr)'bsp_env2.sub_rhs + dual.(bsp_env2.cconstr)'master_env.model[:x]))
            # else
            #     g₂ = Inf
            #     @error "Wrong status2 = $status2"
            # end
        else
            g₂ = 0
        end
        # _UB2 = min(_UB2, v̂₀*g₂ - v̂ₜ)
        _UB2 = g₂ 


        ##################### LB and UB #####################
        LB = τ̂
        # UB = update_UB!(UB,_sx,k̂₀*g₁,v̂₀*g₂,t̂, pConeType)


        @info "Iteration $k: LB = $LB, UB = $UB, _UB1 = $_UB1, _UB2 = $_UB2"


        ##################### check termination #####################
        if  (1e-3 >= _UB1 && 1e-3 >= _UB2 ) || k >= 30 #((UB - LB)/abs(UB) <= 1e-3 || (1e-3 >= _UB1 && 1e-3 >= _UB2 )) || (UB - LB) <= 0.01 ||
            main_env.ifsolved = true
            break
        end

        if time_limit <= time() - start_time
            main_env.ifsolved = true
            break
        end

        ##################### add cuts into DCGLP and master problem #####################
        if termination_status(bsp_env1.model) !=  TIME_LIMIT
            if k̂₀ == 0
                # if status2 == FEASIBLE_POINT
                #     @constraint(main_env.model, main_env.model[:kₜ] >= conπpoints2[end][1]*main_env.model[:k₀] + conπpoints2[end][2]'main_env.model[:kₓ])
                # elseif status2 == INFEASIBILITY_CERTIFICATE
                #     @constraint(main_env.model, 0 >= conπrays2[end][1]*main_env.model[:k₀] + conπrays2[end][2]'main_env.model[:kₓ])
                # end
            else
                # if status1 == FEASIBLE_POINT 
                    # if 1e-3 < _UB1
                    # @info (conπpoints1[end][1] - conπpoints1[end][2])*main_env.model[:k₀]
                    # @info conπpoints1[end][3]'main_env.model[:kₓ]
                    # @info conπpoints1[end][4]*main_env.model[:kₜ]
                        @constraint(main_env.model, 0>= sum(conπpoints1[end][1])*main_env.model[:k₀] - sum(conπpoints1[end][2])*main_env.model[:k₀] + conπpoints1[end][3]'main_env.model[:kₓ] - conπpoints1[end][4]*main_env.model[:kₜ])
                        # @constraint(main_env.model, main_env.model[:vₜ] >= conπpoints1[end][1]*main_env.model[:v₀] + conπpoints1[end][2]'main_env.model[:vₓ])
                        @info "add feasible cut 1"
                    # end
                    # push!(masterconπpoints1, @expression(master_env.model, - master_env.model[:t] + conπpoints1[end][1] + conπpoints1[end][2]'master_env.model[:x]))
                # elseif status1 == INFEASIBILITY_CERTIFICATE
                #     @constraint(main_env.model, 0 >= conπrays1[end][1]*main_env.model[:k₀] + conπrays1[end][2]'main_env.model[:kₓ])
                #     # @constraint(main_env.model, 0 >= conπrays1[end][1]*main_env.model[:v₀] + conπrays1[end][2]'main_env.model[:vₓ])
                #     # push!(masterconπrays1, @expression(master_env.model, conπrays1[end][1] + conπrays1[end][2]'master_env.model[:x]))
                #     @info "add infeasible cut 1"
                # end
            end
        end

        if termination_status(bsp_env2.model) !=  TIME_LIMIT
            if v̂₀ == 0
                # if status1 == FEASIBLE_POINT
                #     @constraint(main_env.model, main_env.model[:vₜ] >= conπpoints1[end][1]*main_env.model[:v₀] + conπpoints1[end][2]'main_env.model[:vₓ])
                # elseif status1 == INFEASIBILITY_CERTIFICATE
                #     @constraint(main_env.model, 0 >= conπrays1[end][1]*main_env.model[:v₀] + conπrays1[end][2]'main_env.model[:vₓ])
                # end
            else
                # if status2 == FEASIBLE_POINT
                    # if 1e-3 < _UB2
                        @constraint(main_env.model, 0 >= sum(conπpoints2[end][1])*main_env.model[:v₀] - sum(conπpoints2[end][2])*main_env.model[:v₀] + conπpoints2[end][3]'main_env.model[:vₓ] - conπpoints2[end][4]*main_env.model[:vₜ])
                        # @constraint(main_env.model, main_env.model[:kₜ] >= conπpoints2[end][1]*main_env.model[:k₀] + conπpoints2[end][2]'main_env.model[:kₓ])
                        @info "add feasible cut 2"
                    # end
                    # push!(masterconπpoints2, @expression(master_env.model, -master_env.model[:t] + conπpoints2[end][1] + conπpoints2[end][2]'master_env.model[:x]))
                # elseif status2 == INFEASIBILITY_CERTIFICATE
                #     @constraint(main_env.model, 0 >= conπrays2[end][1]*main_env.model[:v₀] + conπrays2[end][2]'main_env.model[:vₓ])
                #     # @constraint(main_env.model, 0 >= conπrays2[end][1]*main_env.model[:k₀] + conπrays2[end][2]'main_env.model[:kₓ])
                #     # push!(masterconπrays2, @expression(master_env.model, conπrays2[end][1] + conπrays2[end][2]'master_env.model[:x]))
                    # @info "add infeasible cut 2"
                # end
            end
        end
        # if k>=2
        #     delete(main_env.model, main_env.model[:trust_region])
        #     unregister(main_env.model, :trust_region)
        # end
    end

    ##################### update #####################
    main_env.masterconπpoints1 = masterconπpoints1
    main_env.masterconπpoints2 = masterconπpoints2
    # @info "masterconπpoints1 = $(masterconπpoints1)"
    # @info "masterconπpoints2 = $(masterconπpoints2)"
    main_env.conπpoints1 = conπpoints1  # calculate dual value later need
    main_env.conπpoints2 = conπpoints2

    @info "iteration = $k"
    @info "t̂ = $t̂"
    # @info "k̂₀ = $(k̂₀s[end])"
    # @info "v̂₀ = $(v̂₀s[end])"

    # optimize!(main_env.model)
    # k̂₀ = value(main_env.model[:k₀])
    # k̂ₓ = value.(main_env.model[:kₓ])
    # k̂ₜ = value(main_env.model[:kₜ])
    # v̂₀ = value(main_env.model[:v₀])
    # v̂ₓ = value.(main_env.model[:vₓ])
    # v̂ₜ = value(main_env.model[:vₜ])
    # γₜ = dual(main_env.model[:cont])
    # γ₀ = dual(main_env.model[:con0])
    # γₓ = dual.(main_env.model[:conx])
    # value_κ = -γₜ*k̂ₜ/k̂₀ - γₓ'k̂ₓ/k̂₀ - γ₀
    # value_ν = -γₜ*v̂ₜ/v̂₀ - γₓ'v̂ₓ/v̂₀ - γ₀
    # @info "value_κ = $value_κ"
    # @info "value_ν = $value_ν"
end

function update_UB!(UB,_sx,g₁,g₂,t̂,::L1GammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], Inf)) end
function update_UB!(UB,_sx,g₁,g₂,t̂,::L2GammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], 2)) end
function update_UB!(UB,_sx,g₁,g₂,t̂,::LInfGammaNorm) return min(UB,norm([ _sx; g₁+g₂-t̂], 1)) end






