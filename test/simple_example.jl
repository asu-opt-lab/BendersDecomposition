using JuMP, Gurobi, LinearAlgebra
# define the problem
c = 1
d = [1, 1]
A = [1, 1, 3, 2]
B = [2 -1; 0 -1; -2 3; 1 -1]
b = [3, 2, 3, 2]
function master(c)
    m = Model(Gurobi.Optimizer)
    set_optimizer_attribute(m, MOI.Silent(),true)
    num = length(c)
    @variable(m, x >= 0, Int)
    @variable(m, t >= 0)
    @objective(m, Min, c*x+t)
    return m
end

function sub(d, A, B, b)
    n = size(B, 2)
    nx = size(A, 2)
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, MOI.Silent(),true)
    set_optimizer_attribute(model, "InfUnbdInfo", 1)
    @variable(model, y[1:n] >= 0)
    @variable(model, x)
    @constraint(model, con, A.*x + B*y .>= b)
    @constraint(model, conx, x==0)
    @objective(model, Min, d'y)
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
    while true
        optimize!(master_problem)
        x̂ = value.(master_problem[:x])
        t = value(master_problem[:t])
        LB = JuMP.objective_value(master_problem)


        set_normalized_rhs(sub_problem[:conx], x̂)
        optimize!(sub_problem)
        status = dual_status(sub_problem)
        if status == FEASIBLE_POINT
            subObjVal = JuMP.objective_value(sub_problem)
            extreme_point = dual(sub_problem[:conx])
            const_term = subObjVal - extreme_point * x̂
            ex = @expression(master_problem, -master_problem[:t] + subObjVal + extreme_point * (master_problem[:x] - x̂) )
        elseif status == INFEASIBILITY_CERTIFICATE
            if has_duals(sub_problem)
                subObjVal = JuMP.objective_value(sub_problem)
                extreme_ray = dual(sub_problem[:conx])
                const_term = dual.(sub_problem[:con])'b
                ex = @expression(master_problem, extreme_ray * master_problem[:x] + const_term)
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
        if gap < 1e-3
            break
        end

        if status == FEASIBLE_POINT
            push!(extreme_points, [extreme_point, const_term])
        else
            push!(extreme_rays, [extreme_ray, const_term])
        end

        @constraint(master_problem, 0 >= ex)
    end




    return extreme_points, extreme_rays
end
# master_problem = master(c)
# sub_problem = sub(d, A, B, b)
# @info master_problem, sub_problem
# model = Model();
# @variable(model, x[i=1:2])
# A = [1 2; 3 4]
# b = [5, 6]
# @constraint(model, con_vector, A * x .== b)
run_Benders(c, d, A, B, b)
# n = size(B, 2)
# nx = size(A, 2)
# m = Model(Gurobi.Optimizer)
# set_optimizer_attribute(m, MOI.Silent(),true)
# @variable(m, y[1:n] >= 0)
# @variable(m, x)
# @constraint(m, con, A.*x + B*y .>= b)
