```@meta
CurrentModule = BendersX
```

# [Getting Started](@id getting-started)

This page shows the shortest complete BendersX workflow based on a tiny
capacitated facility location instance. It uses [`BendersSeq`](@ref) and
`HiGHS`, so it does **not** depend on solver callbacks.

## Minimal sequential example

```julia
using BendersX
import BendersX: CFLPData
using JuMP, HiGHS

data = CFLPData(
    2,
    3,
    [4.0, 5.0],
    [2.0, 1.0, 2.0],
    [8.0, 7.0],
    [1.0 3.0 2.0;
     2.0 1.0 2.5],
)

function customize_master_model!(model::Model, data::CFLPData)
    set_optimizer(model, HiGHS.Optimizer)
    set_silent(model)

    @variable(model, x[1:data.n_facilities], Bin)
    @variable(model, t >= -1e6)
    @objective(model, Min, data.fixed_costs' * x + t)
    @constraint(
        model,
        capacity,
        sum(data.capacities[i] * x[i] for i in 1:data.n_facilities) >= sum(data.demands),
    )

    return (x = x,), t
end

function customize_sub_model!(model::Model, data::CFLPData, scen_idx::Int; x)
    set_optimizer(model, HiGHS.Optimizer)
    set_silent(model)

    I, J = data.n_facilities, data.n_customers
    @variable(model, y[1:I, 1:J] >= 0)
    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:, j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.demands .* y[i, :]) <= data.capacities[i] * x[i])

    return nothing
end

master = Master(data; customize = customize_master_model!)
oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
env = BendersSeq(master, oracle; param = BendersSeqParam(verbose = false))
log = solve!(env)
```

## What happened

1. [`Master`](@ref) built the first-stage JuMP model and recorded the coupling
   variables `x` plus the approximation variable `t`.
2. [`ClassicalOracle`](@ref) built an LP-compatible subproblem and reused it at
   each candidate master point.
3. [`BendersSeq`](@ref) orchestrated the solve loop and returned a per-iteration
   log as a `DataFrame`.

## Solver notes

- The example above works with `HiGHS` because it uses a sequential
  environment.
- If you switch to [`BendersBnB`](@ref), [`LazyCallback`](@ref), or
  [`UserCallback`](@ref), you need a solver that supports callbacks.
- The built-in problem helpers in `src/problems/` often configure a different
  optimizer internally. Treat them as templates and adjust the optimizer to fit
  your setup.

## Where to go next

- For the full model-building contract, see [Modeling Guide](@ref modeling-guide).
- For a benchmark-backed example using [`read_cflp_benchmark_data`](@ref), see
  [CFLP Demo](@ref cflp-demo).
- For built-in problem readers and artifact-backed datasets, see
  [Problem Library](@ref problem-library).
