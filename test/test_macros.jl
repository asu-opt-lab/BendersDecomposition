using BendersDecomposition
using JuMP
using CPLEX

facility_cost = [100.0, 120.0]
demand = [10.0, 15.0, 20.0]
shipping_cost = [
    8.0  6.0  3.0;
    4.0  7.0  5.0
]

@benders_decomposition env begin
    @master_problem master begin
        set_optimizer(master, CPLEX.Optimizer)
        set_optimizer_attribute(master, MOI.Silent(), true)
        @variable(master, x[1:2], Bin)
        @variable(master, t >= 0)
        @objective(master, Min, sum(facility_cost[i] * x[i] for i in 1:2) + t)
        @constraint(master, sum(x[i] for i in 1:2) >= 1)
    end

    @sub_problem sub begin
        set_optimizer(sub, CPLEX.Optimizer)
        set_optimizer_attribute(sub, MOI.Silent(), true)
        @variable(sub, x[1:2] >= 0)
        @variable(sub, y[1:2, 1:3] >= 0)
        @objective(sub, Min, sum(shipping_cost[i,j] * y[i,j] for i in 1:2, j in 1:3))
        @constraint(sub, [j=1:3], sum(y[i,j] for i in 1:2) >= demand[j])
        @constraint(sub, [i=1:2], sum(y[i,j] for j in 1:3) <= 1000 * x[i])
    end
end

solver = "CPLEX"
loop_strategy = Sequential()
cut_strategy = ClassicalCut()
params = BendersParams(
    600.0,
    0.00001,
    solver,
    Dict("solver" => solver),
    Dict("solver" => solver),
    Dict(),
    true
)
solve!(env, loop_strategy, cut_strategy, params)