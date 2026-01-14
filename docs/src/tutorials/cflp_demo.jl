# # CFLP Demo
#
# This example demonstrates how to solve the Capacitated Facility Location Problem
# using BendersX.jl with sequential Benders decomposition.

using BendersX
using JuMP, HiGHS

# ## Define the Master Model

function customize_master_model!(model::Model, data::CFLPData)
    optimizer = optimizer_with_attributes(
        HiGHS.Optimizer, MOI.Silent() => true)
    set_optimizer(model, optimizer)

    I = data.n_facilities
    @variable(model, x[1:I], Bin)
    @variable(model, t >= -1e6)
    @objective(model, Min, data.fixed_costs' * x + t)
    @constraint(model, capacity, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))

    return (x = x, ), t
end

# ## Define the Subproblem Model

function customize_sub_model!(model::Model, data::CFLPData, scen_idx::Int; x)
    optimizer = optimizer_with_attributes(
        HiGHS.Optimizer, MOI.Silent() => true)
    set_optimizer(model, optimizer)

    I, J = data.n_facilities, data.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.demands[:] .* y[i, :]) <= data.capacities[i] * x[i])
end

# ## Solve the Problem

data   = read_cflp_benchmark_data("p1")
master = Master(data; customize = customize_master_model!)
oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
env    = BendersSeq(master, oracle)
log    = solve!(env)
