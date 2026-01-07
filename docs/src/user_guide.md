# [User Guide](@id user-guide)

## Design Philosophy
BendersX.jl is built on the principle that Benders decomposition consists of **clear, separable components**, each with a distinct role. By making these components explicit, the framework enables **composable algorithm design**.

The letter "**X**" intentionally has no fixed meaning: it represents *any* extension beyond classical Benders decomposition, including but not limited to:
- alternative Benders cuts,
- stabilization and acceleration techniques,
- different execution and coordination mechanisms, and
- new application domains

## Components
BendersX.jl decomposes the Benders decomposition algorithm into three **explicit and separable components**: the **Master**, the **Oracle**, and the **Environment**.
Each component has a well-defined responsibility and can be independently replaced or extended.

### Master
The **Master** represents the master problem in a Benders decomposition and is responsible for proposing candidate solutions.

Conceptually, the Master:
- represents the current relaxation of the Benders reformulation
- refines this relaxation by incorporating newly generated Benders cuts, and
- produces candidate solutions for evaluation by the Oracle.

In BendersX.jl, the Master is:
- is encapsulated by a subtype of `AbstractMaster`
- owns a JuMP model defining the master problem, and
- accepts newly generated Benders cuts during the solution process.

Users specify the master formulation via a customization hook:
```julia
master = Master(data; customize = customize_master_model!)
```
The function `customize_master_model!` is user-defined where variables, constraints, objectives, and the solver can be specified using standard JuMP syntax. This design ensures that **no modeling language** is introduced beyond JuMP, preserving familiarity and flexibility.

See [Modeling Interface](@ref modeling-interface) for details.

### Oracle
An **Oracle** encapsulates all procedures related to **cut generation** at a given separation point.

In BendersX.jl, an Oracle:
- is implemented as a subtype of `AbstractOracle`
- is fully decoupled from the Master and Environment logic, and
- focuses exclusively on cut generation.

Built-in oracle types include:
- `ClassicalOracle`
- `UnifiedOracle`
- `ParetoOracle`
- `KnapsackOracle` (for facility location)
- `SplitOracle`

Users customize the subproblem formulation via a modeling hook:
```julia
oracle = XOracle(data; customize = customize_sub_model!)
```
The function `customize_sub_model!` defines the JuMP model for the subproblem(s) and is responsible for specifying variables, constraints, objectives, and solver options.
See the [Modeling Interface](@ref modeling-interface) for examples.

Oracle behavior can also be tuned via parameters, such as violation tolerances.

Advanced users can implement custom cut generators by defining a new oracle:
```julia
struct MyOracle <: AbstractOracle
    # fields
end
```
together with the required interface:
```julia
function generate_cuts(oracle::MyOracle)::Tuple
    # return a collection of cuts
end
```
This interface-based design allows new theoretical ideas to be prototyped directly within the framework.

### Environment
The **Environment** controls the **execution logic** of the Benders algorithm.

While the Master and Oracle define what problems are solved and how cuts are generated, the Environment defines *how the algorithm proceeds*. Specifically, it:
- orchestrates the interaction between the Master and Oracle,
- manages iteration order and termination criteria, and
- handles logging, statistics, and output.
In BendersX.jl, an Environment:
- is implemented as a subtype of `AbstractBendersEnv`
- encapsulates the overall control flow of the algorithm, and
- enables alternative execution strategies to be expressed cleanly.

Built-in environments include:
- `BendersSeq` — classical sequential Benders
- `BendersSeqInOut` — sequential Benders with in–out stabilization
- `BendersBnB` — branch-and-bound with Benders cuts

Environment behavior can be adjusted via parameters (e.g., stopping rules, stabilization dynamics).

Custom execution strategies can be implemented by defining a new environment:
```julia
struct MyEnv <: AbstractBendersEnv
    # fields
end
```
together with
```julia
function solve!(env::MyEnv)::DataFrame
    # execution logic
end
```

## Hierarchical Structure of Components
Each component follows a clear type hierarchy that supports specialization and reuse.
This hierarchy enables advanced users to extend existing implementations incrementally rather than starting from scratch.
![dd](OracleHierarchy.pdf)
![dd](EnvHierarchy.pdf)

