# # Stochastic Capacitated Facility Location Problem

# ## Formulation

# Consider the SAA problem:

# ```math
# \begin{aligned}
# \min \ & \sum_{i \in [I]}f_{i} x_{i} + \frac{1}{K} \sum_{k \in [K]} \sum_{i \in [I]}\sum_{j \in [J]} c_{ij}d_{i}y_{ijk} \\
# \text{s.t.} \ & \sum_{i \in [I]} y_{ijk} \geq 1, \ \forall j \in [J], k \in [K] \\
# & \sum_{j \in [J]} y_{ijk} \leq u_{i}x_i, \ \forall i \in [I], k \in [K] \\
# & y_{ijk} \leq x_i, \ \forall i \in [I], j \in [J], k \in [K] \\
# & x \in \mathbb B^{I}, y \in \mathbb{R}_{+}^{I \times J \times K}
# \end{aligned}
# ```

# ## Example

using BendersDecomposition
using CPLEX
using JuMP

# Let's solve this simple instance directly using a solver.

data = read_stochastic_capacited_facility_location_problem("f10-c10-s25-r3-1"; filepath=joinpath(dirname(dirname(dirname(@__DIR__))), "data", "SCFLP"))
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
cut_strategy = ClassicalCut()
solution_procedure = StochasticSequential()
run_Benders(data, solution_procedure, cut_strategy, params)
#-

# If you want to use knapsack cuts, you can use the following code:
cut_strategy = KnapsackCut()
solution_procedure = StochasticSequential()
run_Benders(data, solution_procedure, cut_strategy, params)
#-

# If you want to use callback methods, you can use the following code:
cut_strategy = ClassicalCut()
solution_procedure = StochasticCallback()
run_Benders(data, solution_procedure, cut_strategy, params)
#-