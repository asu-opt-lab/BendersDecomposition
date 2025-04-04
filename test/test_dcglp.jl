using Test
using BendersOracle
using CPLEX
using Gurobi
using JuMP


facility_cost = [100.0, 120.0]
demand = [10.0, 15.0, 20.0]
shipping_cost = [
    8.0  6.0  3.0;
    4.0  7.0  5.0
]

# @disjunctive_benders_decomposition env begin
#     # @coupling_variable(env, x[1:2], Bin)
#     @master_problem master begin
#         # set_optimizer(master, CPLEX.Optimizer)
#         # set_optimizer_attribute(master, MOI.Silent(), true)
#         @variable(master, x[1:2] >= 0, Bin)
#         @variable(master, t >= 0)
#         @objective(master, Min, sum(facility_cost[i] * x[i] for i in 1:2) + t)
#         @constraint(master, sum(x[i] for i in 1:2) >= 1)
#     end

#     @sub_problem sub begin
#         # set_optimizer(sub, CPLEX.Optimizer)
#         # set_optimizer_attribute(sub, MOI.Silent(), true)
#         @variable(sub, x[1:2])
#         @variable(sub, y[1:2, 1:3] >= 0)
#         @objective(sub, Min, sum(shipping_cost[i,j] * y[i,j] for i in 1:2, j in 1:3))
#         @constraint(sub, [j=1:3], sum(y[i,j] for i in 1:2) >= demand[j])
#         @constraint(sub, [i=1:2], sum(y[i,j] for j in 1:3) <= 1000 * x[i])
#     end
# end

@disjunctive_benders_decomposition env begin
    # @coupling_variable(env, x[1:2], Bin)
    @master_problem master begin
        set_optimizer(master, CPLEX.Optimizer)
        set_optimizer_attribute(master, MOI.Silent(), true)
        @variable(master, x[1:2] >= 0, Bin)
        @variable(master, t >= 0)
        @objective(master, Min, sum(facility_cost[i] * x[i] for i in 1:2) + t)
        @constraint(master, sum(x[i] for i in 1:2) >= 1)
    end

    @sub_problem sub begin
        set_optimizer(sub, Gurobi.Optimizer)
        set_optimizer_attribute(sub, MOI.Silent(), true)
        set_optimizer_attribute(sub, "InfUnbdInfo", 1)
        @variable(sub, x[1:2])
        @variable(sub, y[1:2, 1:3] >= 0)
        @objective(sub, Min, sum(shipping_cost[i,j] * y[i,j] for i in 1:2, j in 1:3))
        @constraint(sub, [j=1:3], sum(y[i,j] for i in 1:2) >= demand[j])
        @constraint(sub, [i=1:2], sum(y[i,j] for j in 1:3) <= 1000 * x[i])
    end

    @config DisjunctiveCut(
        ClassicalCut(),
        L1Norm(),
        PureDisjunctiveCut(),
        true,
        false,
        false,
        true
    )
end 
# println(env.master)

# println(env.sub.model)

# println(env.sub.fixed_x_constraints)

println(env.sub.other_constraints)

# println(env.sub.fixed_x_values)

# println(env.dcglp.model)

set_optimizer(env.dcglp.model, CPLEX.Optimizer)
set_optimizer_attribute(env.dcglp.model, MOI.Silent(), true)
loop_strategy = GenericSequential(10.0, 5, 1e-6, true)

cut_strategy = DisjunctiveCut(
    ClassicalCut(),
    L1Norm(),
    PureDisjunctiveCut(),
    true,
    false,
    false,
    true
)

solve!(env, loop_strategy, cut_strategy)
