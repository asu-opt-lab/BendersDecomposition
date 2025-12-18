*End-to-end workflow for solving the Capacitated Facility Location Problem (CFLP) using built-in BendersX.jl components*

## 1. Mathematical Modeling
### [Parameters](@id cflp-parameter)
- number of facilities ``I``
- number of customers ``J``
- facility capacities ``\mathbf{b} = (b_i)_{i \in [I]}``
- customer demands ``\mathbf{d} = (d_j)_{j \in [J]}``
- fixed opening cost ``\mathbf{f} = (f_i)_{i \in [I]}``
- transportation cost ``\mathbf{c} = (c_{ij})_{i \in [I], j \in [J]}``. 

The CFLP can be written as the following mixed-integer program:
```math 
    \begin{align}
    \min \ & \sum_{i \in [I]}f_ix_i + \sum_{i \in [I], j \in [J]}c_{ij}d_j y_{ij}\nonumber \\
    \text{s.t.} \ & \sum_{i \in [I]} y_{ij } = 1, \ \forall j \in [J],\nonumber \\
    & y_{ij} \le x_i, \ \forall i \in [I], j \in [J],\nonumber \\
    & \sum_{j \in [J]}d_j y_{ij} \le b_i x_i, \ \forall i \in [I]\nonumber \\
    & y_{ij} \ge 0, \ \forall i \in [I], j \in [J],\nonumber \\
    &x_i \in \mathbb B, \ \forall i \in [I],\nonumber 
\end{align}
```
- The first constraint enforces demand fulfillment.
- The second links assignment decisions to facility-opening decisions for better LP relaxation.
- The third enforces facility capacity limits.

### Benders Reformulation 
Separating the binary variables and introducing an auxiliary variable ``t``, we obtain
#### [Master problem](@id cflp-master)
```math 
    \begin{align}
    \min \ & \mathbf{f}^\top \mathbf{x} + t \nonumber \\
    \text{s.t.} \ 
    &x_i \in \mathbb B, \ \forall i \in [I]. \nonumber 
    \end{align}
```
Two redundant constraints may optionally be added:
```math 
\begin{align}
     & \mathbf{b}^\top \mathbf{x} \ge \mathbf{1}^\top \mathbf{d}, \nonumber \\
    & t \ge -10^{6}, \nonumber
\end{align}
```
where the first ensures sufficient total capacity and the second prevents 
``t`` from initially taking an arbitrarily negative value.

#### [Subproblem](@id cflp-sub)
```math
\begin{align}
    \min \ & \sum_{i \in [I], j \in [J]}c_{ij}d_j y_{ij} \nonumber \\
    \text{s.t.} \ & \sum_{i \in [I]} y_{ij } = 1, \ \forall j \in [J],\nonumber \\
    & y_{ij} \le x_i, \ \forall i \in [I], j \in [J],\nonumber \\
    & \sum_{j \in [J]}d_j y_{ij} \le b_i x_i, \ \forall i \in [I]\nonumber \\
    & y_{ij} \ge 0, \ \forall i \in [I], j \in [J].\nonumber \\
\end{align}
```




## 2. Full Code Example
The following code block illustrates how to solve a CFLP instance using a sequential Benders decomposition with classical optimality and feasibility cuts.

```julia
using BendersX
using JuMP, CPLEX

struct CFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    capacities::Vector{Float64}
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end

function read_cflp_benchmark_data(path_to_raw_data)
    """
    Users are responsible for loading their raw data into the 
    user-defined data structure (e.g., `CFLPData`). See, for example,
    `read_cflp_benchmark_data` in
    `./src/BendersLibrary/src/problems/cflp/data_reader.jl`.
    """
end

function customize_master_model!(model::Model, data::CFLPData)
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, MOI.Silent() => true)
    set_optimizer(model, optimizer)

    I = data.n_facilities
    @variable(model, x[1:I], Bin)
    @variable(model, t >= -1e6)
    @objective(model, Min, data.fixed_costs'* x + t)
    @constraint(model, capacity, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))

    return (x = x, ), t
end

function customize_sub_model!(model::Model, data::CFLPData, scen_idx::Int; x) 
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, MOI.Silent() => true)
    set_optimizer(model, optimizer)

    I, J = data.n_facilities, data.n_customers   
    @variable(model, y[1:I, 1:J] >= 0)
    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.demands[:] .* y[i,:]) <= data.capacities[i] * x[i])
end

data   = read_cflp_benchmark_data("p1")
master = Master(data; customize = customize_master_model!)
oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
env    = BendersSeq(master, oracle)
log    = solve!(env)
```




## 2. Data 
*Providing Instance Data to the Benders Engine*

Users must define a subtype of `AbstractData` to store all problem-specific parameters required by the master and subproblem models. `BendersX.jl` does not impose any restrictions on the structure or fields of this type—any user-defined container is acceptable as long as it provides the information needed to build the models.

### Example
The [parameters](@ref cflp-parameter) of the Capacitated Facility Location Problem (CFLP) may be written as:
```julia
struct CFLPData <: AbstractData
    n_facilities::Int
    n_customers::Int
    capacities::Vector{Float64}
    demands::Vector{Float64}
    fixed_costs::Vector{Float64}
    costs::Matrix{Float64}
end
```
### Loading Data
Users are responsible for loading and preprocessing raw instance data into their `AbstractData` subtype. A typical pattern is:
```julia
read_data(path_to_raw_data) -> MyData
```
For an example specific to CFLP, see `read_cflp_benchmark_data` in `./src/BendersLibrary/src/problems/cflp/data_reader.jl`.

## [3. Modeling Interface](@id modeling-interface)
*Defining Master and Subproblem Models in BendersX.jl*

Users provide the master and subproblem formulations through *customization functions* written in standard JuMP syntax.
If you are unfamiliar with JuMP, please refer to the JuMP.jl documentation for an introduction:
[Julia JuMP](https://jump.dev/JuMP.jl/stable/)

### Master Modeling
Users specify the master formulation by implementing a function of the form:
```julia
customize_master_model!(model::Model, data::AbstractData) -> NamedTuple, Vector{VariableRef}
```
Within this function, users use standard JuMP commands to declare master-level variables, constraints, and the objective, as well as configure the optimizer.
The function must return:
1. a NamedTuple mapping symbolic variable names to non-auxiliary master variables, and
2. a Vector{VariableRef} containing the auxiliary variables $t$ used for Benders cuts.

All modeling and solver-related decisions are entirely under the user’s control, while the Benders engine itself remains independent of the model and the solver.

### Example
[The CFLP master problem](@ref cflp-master) can be implemented as:

```julia
function customize_master_model!(model::Model, data::CFLPData)
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, MOI.Silent() => true)
    set_optimizer(model, optimizer)

    I = data.n_facilities
    @variable(model, x[1:I], Bin)
    @variable(model, t >= -1e6)
    @objective(model, Min, data.fixed_costs'* x + t)
    @constraint(model, capacity, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))

    return (x = x, ), t
end
```

### Subproblem Modeling
Subproblems are specified by the user through:
```julia
customize_sub_model!(model::Model, data::AbstractData, scen_idx::Int; kwargs...)
```
Here, `kwargs...`` contains the symbolic names of the master variables that appear in the subproblem. This allows users to formulate the subproblem in JuMP **while referencing these master variables directly**, without explicitly adding them to the subproblem model.

### Example
[The CFLP subproblem](@ref cflp-sub) can be implemented like this:
```julia
function customize_sub_model!(model::Model, data::CFLPData, scen_idx::Int; x) 
    optimizer = optimizer_with_attributes(
        CPLEX.Optimizer, "CPXPARAM_Threads" => 7, MOI.Silent() => true)
    set_optimizer(model, optimizer)

    I, J = data.n_facilities, data.n_customers   
    @variable(model, y[1:I, 1:J] >= 0)
    cost_demands = data.costs .* data.demands'
    @objective(model, Min, sum(cost_demands .* y))
    @constraint(model, demand[j in 1:J], sum(y[:,j]) == 1)
    @constraint(model, facility_open, y .<= x)
    @constraint(model, capacity[i in 1:I], sum(data.demands[:] .* y[i,:]) <= data.capacities[i] * x[i])
end
```

## 4. Running Benders

Users can employ suitable subtypes of `AbstractOracle` (e.g., `UnifiedOracle`, `ParetoOracle`, `CFLKnapsackOracle` for the CFLP) and any subtype of `AbstractBendersEnv` provided by BendersX.jl (see Figures \ref{fig:env} and \ref{fig:oracle}).

Furthermore, users can configure both the environment and the oracle through dedicated parameter objects or customizable sub-components.

How to configure environment and oracles. 


