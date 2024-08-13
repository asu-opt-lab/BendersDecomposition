using JuMP, Gurobi, LinearAlgebra, Plots
# using PlotlyJS
# define the problem
c = [1,1]
d = 1
A = [15 10 10 -10 -70; 45 7 16 35 26]
# B = [1, 3, 1, 2, 2]
# B = [2, 3, 4, 2, 5] # correct 
B = [10, 6, 4, 7, 9]
b = [8, 17, 9, 1, 49]

# function mip(d, A, B, b, c)
#     m = Model(Gurobi.Optimizer)
#     set_optimizer_attribute(m, MOI.Silent(),true)
#     @variable(m, x, Bin)
#     @variable
#     return m
# end
function master(c)
    m = Model(Gurobi.Optimizer)
    set_optimizer_attribute(m, MOI.Silent(),true)
    num = length(c)
    # @variable(m, x, Bin)
    @variable(m, 0<=x[1:2]<=1)
    @variable(m, t>=-1e06)
    @objective(m, Min, c'x+t)
    # @constraint(m, t +15x >= 8)
    # @constraint(m, t -35x >= -24.5)
    return m
end

function sub(d, A, B, b)
    n = size(B, 2)
    nx = size(A, 2)
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, MOI.Silent(),true)
    # set_optimizer_attribute(model, "InfUnbdInfo", 1)
    @variable(model, y >= 0)
    @variable(model, x[1:2])
    @constraint(model, con, A'x + B.*y .>= b)
    @constraint(model, conx, x.==0)
    @objective(model, Min, d*y)
    return model
end

function _DCGLP()

    # model = Model(CPLEX.Optimizer)
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, MOI.Silent(),true)

    
    # Variables
    @variable(model, τ)
    @variable(model, k₀>=0)
    @variable(model, kₓ[1:2])
    @variable(model, kₜ)
    @variable(model, v₀>=0)
    @variable(model, vₓ[1:2])
    @variable(model, vₜ)
    @variable(model, sx)
    @variable(model, st)

    # Objective
    @objective(model, Min, τ)
    

    # Constraints
    @constraint(model, consigma1, 0 .>= 1*k₀ .- kₓ) 
    @constraint(model, coneta1, 0 .>= -k₀ .+ kₓ) 
    @constraint(model, consigma2, 0 .>= -0*v₀ .+ vₓ) 
    @constraint(model, coneta2, 0 .>= -v₀ .+ vₓ)

    @constraint(model, conv1, 0 .>= -kₓ)
    @constraint(model, conv2, 0 .>= -vₓ)

    @constraint(model, con0, k₀ + v₀ == 1)
    
    # @constraint(model, concone, [τ; sx ; st] in MOI.NormInfinityCone(2 + 1))
    @constraint(model, concone, [τ; sx ; st] in MOI.NormOneCone(2 + 1))

    
    @constraint(model, conx, kₓ .+ vₓ .- sx .== 0)  #x̂
    @constraint(model, cont, kₜ + vₜ - st == 0) #t̂
    
    
    # for i in eachindex(sub_env.split_info.γ₀s)
    #     @constraint(model, sub_env.split_info.γ₀s[i]*k₀ + sub_env.split_info.γₓs[i]'kₓ + sub_env.split_info.γₜs[i]*kₜ <= 0)
    #     @constraint(model, sub_env.split_info.γ₀s[i]*v₀ + sub_env.split_info.γₓs[i]'vₓ + sub_env.split_info.γₜs[i]*vₜ <= 0)
    # end

    return model
end

function run_Benders(c, d, A, B, b)
    master_problem = master(c)
    sub_problem = sub(d, A, B, b)
    num = length(c)
    t = 0
    LB = -Inf
    UB = Inf
    extreme_points = []
    extreme_rays = []
    iter = 1
    while true
        optimize!(master_problem)
        x̂ = value.(master_problem[:x])
        t = value(master_problem[:t])
        LB = JuMP.objective_value(master_problem)
        @info "x̂ = $x̂"

        set_normalized_rhs.(sub_problem[:conx], x̂)
        optimize!(sub_problem)
        @info value(sub_problem[:y])
        status = dual_status(sub_problem)
        if status == FEASIBLE_POINT
            subObjVal = JuMP.objective_value(sub_problem)
            extreme_point = dual.(sub_problem[:conx])
            const_term = subObjVal - extreme_point'x̂
            ex = @expression(master_problem, -master_problem[:t] + subObjVal + extreme_point'(master_problem[:x] .- x̂) )
        elseif status == INFEASIBILITY_CERTIFICATE
            if has_duals(sub_problem)
                subObjVal = JuMP.objective_value(sub_problem)
                extreme_ray = dual.(sub_problem[:conx])
                const_term = dual.(sub_problem[:con])'b
                ex = @expression(master_problem, extreme_ray'master_problem[:x] + const_term)
            else
                @error "infeasible sub has no infeasible ray"
                throw(-1)
            end      
            subObjVal = 1e+99
        else
            @error "dual of sub is neither feasible nor infeasible certificate: $status"
            # throw(-1)
            ex = 0
            subObjVal = -Inf
        end

        UB = min(UB, c*x̂ + subObjVal)
        gap = 100 * (UB - LB)/ abs(UB)
        @info "LB: $LB, UB: $UB, Gap: $gap"
        if gap < 1e-3 || iter >= 5
            break
        end

        if status == FEASIBLE_POINT
            push!(extreme_points, [extreme_point, const_term])
        else
            push!(extreme_rays, [extreme_ray, const_term])
        end

        @constraint(master_problem, 0 >= ex)
        iter += 1
    end



    @info objective_value(master_problem)
    return extreme_points, extreme_rays
end

function txfigure(extreme_points, extreme_rays)
    x = range(0,1, length=10)
    figure = plot()
    for i in 1:length(extreme_points)
        @info i
        extreme_point = extreme_points[i][1]
        const_term = extreme_points[i][2]
        y = const_term .+ extreme_point * x
        plot!(figure, x, y, label="extreme_point_$i")
    end

    savefig("extreme_points.png")
end

function run_Benders_split(c, d, A, B, b)

    figure = plot()
    x1 = range(0,3, length=10)
    x2 = range(0,3, length=10)
    # layout = Layout(
    #     title="Mt Bruno Elevation",
    #     autosize=false,
    #     width=500,
    #     height=500,
    #     margin=attr(l=65, r=50, b=65, t=90)
    # )
    # layout = Layout(
    #     scene=attr(
    #         xaxis=attr(
    #             nticks=4,
    #             range=[0,1]
    #         ),
    #         yaxis=attr(
    #             nticks=4,
    #             range=[-0,1]
    #         ),
    #         zaxis=attr(
    #             nticks=4,
    #             range=[-50,50]
    #         ),
    #     ),
    #     width=700,
    #     margin=attr(
    #         r=20,
    #         l=10,
    #         b=10,
    #         t=10
    #     ),
    # )
    master_problem = master(c)
    sub_problem = sub(d, A, B, b)
    num = length(c)
    t = 0
    LB = -Inf
    UB = Inf
    extreme_points = []
    extreme_rays = []
    iter = 1
    DCGLP = _DCGLP()

    while true
        optimize!(master_problem)
        x̂ = value.(master_problem[:x])
        t̂ = value(master_problem[:t])
        LB = JuMP.objective_value(master_problem)
        @info "x̂ = $x̂"
        @info "master_problem"
        @info master_problem

        set_normalized_rhs.(DCGLP[:conx], x̂)
        set_normalized_rhs.(DCGLP[:cont], t̂)
        _UB1 = Inf  
        _UB2 = Inf
        k = 1
        while true
            @info DCGLP
            optimize!(DCGLP)
            k̂₀ = value(DCGLP[:k₀])
            k̂ₓ = value.(DCGLP[:kₓ])
            k̂ₜ = value(DCGLP[:kₜ])
            v̂₀ = value(DCGLP[:v₀])
            v̂ₓ = value.(DCGLP[:vₓ])
            v̂ₜ = value(DCGLP[:vₜ])
            τ̂ = value(DCGLP[:τ])
            _sx = value(DCGLP[:sx])
            @info "k̂₀ = $k̂₀, v̂₀ = $v̂₀, k̂ₓ = $k̂ₓ, v̂ₓ = $v̂ₓ, k̂ₜ = $k̂ₜ, v̂ₜ = $v̂ₜ"
            if k̂₀ != 0
                set_normalized_rhs.(sub_problem[:conx], k̂ₓ./k̂₀)
                optimize!(sub_problem)
                status1 = dual_status(sub_problem)
                if status1 == FEASIBLE_POINT
                    g₁ = JuMP.objective_value(sub_problem)
                    extreme_point = dual.(sub_problem[:conx])
                    const_term = dual.(sub_problem[:con])'b
                    ex1 = @expression(DCGLP, -DCGLP[:kₜ] + const_term*DCGLP[:k₀] + extreme_point'DCGLP[:kₓ] )
                    # push!(extreme_points, [extreme_point, const_term])
                    _UB1 = g₁ - k̂ₜ
                elseif status1 == INFEASIBILITY_CERTIFICATE
                    g₁ = Inf
                    extreme_ray = dual.(sub_problem[:conx])
                    const_term = dual.(sub_problem[:con])'b
                    ex1 = @expression(DCGLP, const_term*DCGLP[:k₀] + extreme_ray'DCGLP[:kₓ])
                    # push!(extreme_rays, [extreme_ray, const_term])
                else
                    g₁ = Inf
                    @error "Wrong status1 = $status1"
                end
            else
                g₁ = 0
            end
            _UB1 = min(_UB1, k̂₀*g₁ - k̂ₜ)

            if v̂₀ != 0
                set_normalized_rhs.(sub_problem[:conx], v̂ₓ./v̂₀)
                optimize!(sub_problem)
                status2 = dual_status(sub_problem)
                if status2 == FEASIBLE_POINT
                    g₂ = JuMP.objective_value(sub_problem)
                    extreme_point2 = dual.(sub_problem[:conx])
                    const_term2 = dual.(sub_problem[:con])'b
                    ex2 = @expression(DCGLP, -DCGLP[:vₜ] + const_term2*DCGLP[:v₀] + extreme_point2'DCGLP[:vₓ] )
                    _UB2 = g₂ - v̂ₜ
                    # push!(extreme_points, [extreme_point2, const_term2])
                elseif status2 == INFEASIBILITY_CERTIFICATE
                    g₂ = Inf
                    extreme_ray2 = dual.(sub_problem[:conx])
                    const_term2 = dual.(sub_problem[:con])'b
                    ex2 = @expression(DCGLP, const_term2*DCGLP[:v₀] + extreme_ray2'DCGLP[:vₓ])
                    # push!(extreme_rays, [extreme_ray2, const_term2])
                else
                    g₂ = Inf
                    @error "Wrong status2 = $status2"
                end
            else
                g₂ = 0
            end
            _UB2 = min(_UB2, v̂₀*g₂ - v̂ₜ)

            LB = τ̂
            # UB = min(UB,norm([ _sx; g₁+g₂-t̂], Inf))
            UB = min(UB,norm([ _sx; g₁+g₂-t̂], 1))

            @info "Iteration $k: LB = $LB, UB = $UB, _UB1 = $_UB1, _UB2 = $_UB2"
            if ((UB - LB)/abs(UB) <= 1e-6 || (1e-3 >= _UB1 && 1e-3 >= _UB2 )) || (UB - LB) <= 0.01 || k >= 5
                @info DCGLP
                break
            end

            if k̂₀ == 0
                if status2 == FEASIBLE_POINT
                    @constraint(DCGLP, 0 >= -DCGLP[:kₜ] + const_term2*DCGLP[:k₀] + extreme_point2'DCGLP[:kₓ])
                else
                    @constraint(DCGLP, 0 >= const_term2*DCGLP[:k₀] + extreme_ray2'DCGLP[:kₓ])
                end
            else
                # @constraint(DCGLP, 0 >= ex1)
                if status1 == FEASIBLE_POINT
                    if 1e-3 < _UB1
                        @constraint(DCGLP, 0 >= ex1)
                        push!(extreme_points, [extreme_point, const_term])
                        # @constraint(master_problem, master_problem[:t] >= const_term + extreme_point * master_problem[:x])
                    end
                else
                    @constraint(DCGLP, 0 >= ex1)
                    push!(extreme_rays, [extreme_ray, const_term])
                    # @constraint(master_problem, 0 >= const_term + extreme_ray * master_problem[:x])
                end
            end

            if v̂₀ == 0
                if status1 == FEASIBLE_POINT
                    @constraint(DCGLP, 0 >= -DCGLP[:vₜ] + const_term*DCGLP[:v₀] + extreme_point'DCGLP[:vₓ])
                else
                    @constraint(DCGLP, 0 >= const_term*DCGLP[:v₀] + extreme_ray'DCGLP[:vₓ])
                end
            else
                # @constraint(DCGLP, 0 >= ex2)
                if status2 == FEASIBLE_POINT
                    if 1e-3 < _UB2
                        @constraint(DCGLP, 0 >= ex2)
                        push!(extreme_points, [extreme_point2, const_term2])
                        # @constraint(master_problem, master_problem[:t] >= const_term2 + extreme_point2 * master_problem[:x])
                    end
                else
                    @constraint(DCGLP, 0 >= ex2)
                    push!(extreme_rays, [extreme_ray2, const_term2])
                    # @constraint(master_problem, 0 >= const_term2 + extreme_ray2 * master_problem[:x])
                end
            end

            k+=1
        end

        γₜ = dual(DCGLP[:cont])
        γ₀ = dual(DCGLP[:con0])
        γₓ = dual.(DCGLP[:conx])
        @info "γ₀ = $γ₀, γₓ = $γₓ, γₜ = $γₜ"
        @constraint(master_problem, -γ₀ - γₓ'master_problem[:x] - γₜ*master_problem[:t] >= 0) 
        @info master_problem
        f(x1,x2)= -γ₀/γₜ .- γₓ[1]/γₜ*x1 .- γₓ[2]/γₜ*x2
        plot!(figure, x1, x2, f, st=:surface,label="gamma_$iter")
        # @info f
        display(figure)
        # p = plot(surface(zdata=f,x=x1,y=x2), layout)
        # p = plot(mesh3d(x = (x1), y = (x2), z = (f.(x1,x2))), layout)
        # p = plot(mesh3d(
        #     x=(0.7 .* randn(10)),
        #     y=(0.55 .* randn(10)),
        #     z=(0.40 .* randn(10)),
        #     color="rgba(244,22,100,0.6)"
        # ),
        # layout,)
        # p = surface(x=x1, y=x2, z=f)
        # display(p)
        if iter >= 4
            break
        end
        

        # if status == FEASIBLE_POINT
        #     push!(extreme_points, [extreme_point, const_term])
        # else
        #     push!(extreme_rays, [extreme_ray, const_term])
        # end

        # @constraint(master_problem, 0 >= ex)
        iter += 1
    end


    optimize!(master_problem)
    @info objective_value(master_problem)
    for i in 1:length(extreme_points)
        extreme_point = extreme_points[i][1]
        const_term = extreme_points[i][2]
        g(x1,x2) = const_term .+ extreme_point[1] * x1 .+ extreme_point[2] * x2
        plot!(figure, x1, x2, g, st=:surface,label="extreme_point_$i")
    end
    savefig("split_benders_3d.png")
    # figure
    return extreme_points, extreme_rays
end


run_Benders_split(c, d, A, B, b)
# extreme_points, extreme_rays = run_Benders(c, d, A, B, b)
# txfigure(extreme_points, extreme_rays)