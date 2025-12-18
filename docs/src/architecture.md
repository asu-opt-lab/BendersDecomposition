
# Architecture 
BendersX.jl features a modular, hierarchical architecture composed of three main components—**master**, **oracle**, and **environment**—each responsible for a distinct part of the algorithm.

## Master
- Represents the master problem in a Benders decomposition.
- Generates candidate solutions and incorporates newly generated cuts into the model.
- Exposes hooks for user-defined modeling code:
```julia
master = Master(data; customize = customize_master_model!)
```
The master model is specified through the user-provided function `customize_master_model!`; see, for example, [Modeling Interface](@ref modeling-interface). 

### Oracle
![dd](OracleHierarchy.pdf)
- Encapsulates all procedures for cut generation given a separation point.
- Built-in oracle types include:
    - ClassicalOracle
    - UnifiedOracle
    - ParetoOracle
    - KnapsackOracle (for facility location)
    - SplitOracle
- Users can implement customized cut-generation by defining:
```julia
struct MyOracle <: AbstractOracle
    # fields
end
```
together with
```julia
function generate_cuts() -> 
```

### Environment
![dd](EnvHierarchy.pdf)
- Orchestrates the iteration logic of the Benders algorithm.
- Controls how master and oracle interact each iteration.
- Built-in environments include:
    - BendersSeq — classic sequential algorithm
    - BendersSeqInOut — sequential algorithm with in-out stabilization
    - BendersBnB — branch-and-bound with Benders cuts
    - Extensible through:
    ```julia
    struct MyEnv <: AbstractBendersEnv
        # fields
    end
    ```
    together with
    ```julia
    function solve! -> DataFrame
    ```