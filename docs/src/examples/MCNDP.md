```@meta
EditURL = "MCNDP.jl"
```

# Multi-Commodity Network Design Problem

Consider this problem:

```math
\begin{aligned}
[\text{MCNDP}] \quad \min & \sum_{i\in I}\sum_{j\in A}d_ic_{ij}x_{ij} + \sum_{j\in A}f_jy_j \\
\text{s.t.} \quad & \sum_{j\in A_v^+}x_{ij} - \sum_{j\in A_v^-}x_{ij} = b_{iv} & \forall i \in I, v \in V \\
& \sum_{i\in I}d_i x_{ij} \leq s_jy_j & \forall j \in A \\
& x_{ij} \leq y_j & \forall i \in I, j \in A \\
& \mathbf{x} \geq \mathbf{0}, \mathbf{y} \in Y
\end{aligned}
```

## Example

````@example MCNDP
using BendersDecomposition
using Gurobi
using JuMP
````

Let's solve this simple instance directly using a solver.

````@example MCNDP
data = read_mcndp_instance("r01.1.dow"; filepath=joinpath(dirname(dirname(dirname(@__DIR__))), "data", "NDR"))
milp = create_milp(data)
set_optimizer(milp.model, Gurobi.Optimizer)
set_optimizer_attribute(milp.model, MOI.Silent(), false)
optimize!(milp.model)
````

Then, let's solve this problem using Benders decomposition.

````@example MCNDP
solver = "Gurobi"
params = BendersParams(
    60.0, # Time limit
    0.00001, # Tolerance
    solver,
    Dict("solver" => solver),
    Dict("solver" => solver, "InfUnbdInfo" => 1),
    Dict(),
    true # verbose
)
cut_strategy = ClassicalCut()
solution_procedure = Sequential()
run_Benders(data, solution_procedure, cut_strategy, params)
````

If you want to use knapsack cuts, you can use the following code:

````@example MCNDP
cut_strategy = KnapsackCut()
solution_procedure = Sequential()
run_Benders(data, solution_procedure, cut_strategy, params)
````

If you want to use callback methods, you can use the following code:

````@example MCNDP
cut_strategy = ClassicalCut()
solution_procedure = Callback()
run_Benders(data, solution_procedure, cut_strategy, params)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

