```@meta
EditURL = "CFLP.jl"
```

# Capacitated Facility Location Problem

## Formulation
Consider this problem:

```math
\begin{aligned}
\min \ & \sum_{i \in [I]}f_{i} x_{i} + \sum_{i \in [I]}\sum_{j \in [J]} c_{ij}d_{i}y_{ij} \\
\text{s.t.} \ & \sum_{i \in [I]} y_{ij} \geq 1, \ \forall j \in [J] \\
& \sum_{j \in [J]} y_{ij} \leq u_{i}x_i, \ \forall i \in [I] \\
& y_{ij} \leq x_i, \ \forall i \in [I], j \in [J] \\
& x \in \mathbb B^{I}, y \in \mathbb{R}_{+}^{I \times J}
\end{aligned}
```

## Example

````@example CFLP
using BendersDecomposition
using CPLEX
using JuMP
````

Let's solve this simple instance directly using a solver.

````@example CFLP
data = read_cflp_benchmark_data("p1"; filepath=joinpath(dirname(dirname(dirname(@__DIR__))), "data", "locssall"))
milp = create_milp(data)
set_optimizer(milp.model, CPLEX.Optimizer)
set_optimizer_attribute(milp.model, MOI.Silent(), false)
optimize!(milp.model)
````

Then, let's solve this problem using Benders decomposition.

````@example CFLP
solver = "CPLEX"
params = BendersParams(
    60.0, # Time limit
    0.00001, # Tolerance
    solver,
    Dict("solver" => solver),
    Dict("solver" => solver),
    Dict(),
    true # verbose
)
cut_strategy = ClassicalCut()
solution_procedure = Sequential()
run_Benders(data, solution_procedure, cut_strategy, params)
````

If you want to use knapsack cuts, you can use the following code:

````@example CFLP
cut_strategy = KnapsackCut()
solution_procedure = Sequential()
run_Benders(data, solution_procedure, cut_strategy, params)
````

If you want to use callback methods, you can use the following code:

````@example CFLP
cut_strategy = ClassicalCut()
solution_procedure = Callback()
run_Benders(data, solution_procedure, cut_strategy, params)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

