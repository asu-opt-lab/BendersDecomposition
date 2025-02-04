```@meta
EditURL = "SNIP.jl"
```

# Stochastic Network Interdiction Problem

## Formulation
```math
\begin{aligned}
\text{[SNIP]} \quad \min & \sum_{k} p_k x_{i,k} \\
\text{s.t.} \quad & x_{i,k} = 1 & \forall k \in K \\
& x_{i,k} - q_kx_{j,k} \geq 0 & \forall a = (i,j) \in D, k \in K \\
& x_{i,k} - r_kx_{j,k} \geq 0 & \forall a = (i,j) \in A \setminus D, k \in K \\
& x_{i,k} - r_kx_{j,k} \geq -(r_k - q_k)y_ay_k & \forall a = (i,j) \in D, k \in K \\
& x \geq 0, y \in Y
\end{aligned}
```

## Example

````@example SNIP
using BendersDecomposition
using JuMP
using CPLEX
````

# Let's solve this simple instance directly using a solver.

data = read_snip_data(0, 1, 90.0; base_dir=joinpath(dirname(dirname(dirname(@__DIR__))), "data", "SNIP"))
milp = create_milp(data)
set_optimizer(milp.model, CPLEX.Optimizer)
set_optimizer_attribute(milp.model, MOI.Silent(), false)
optimize!(milp.model)

# Then, let's solve this problem using Benders decomposition.
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

# Then, let's solve this problem using Benders decomposition.
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
solution_procedure = StochasticSequential()
run_Benders(data, solution_procedure, cut_strategy, params)
#-

# If you want to use callback methods, you can use the following code:
cut_strategy = ClassicalCut()
solution_procedure = StochasticCallback()
run_Benders(data, solution_procedure, cut_strategy, params)
#-

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

