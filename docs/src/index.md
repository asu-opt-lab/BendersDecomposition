# BendersX.jl Documentation

Welcome to the documentation for **BendersX.jl**.

## Introduction
**BendersX.jl** is a modular and extensible framework for implementing Benders decomposition algorithms in Julia. The package separates the algorithm into three core components—**Master**, **Oracle**, and **Environment**—allowing each part of the algorithm to be customized, replaced, or extended independently. This design enables rapid prototyping, reproducible experiments, and fair comparison of alternative Benders strategies without modifying modeling code.

Built on top of JuMP, BendersX.jl allows users to formulate master and subproblem models using standard JuMP modeling syntax while delegating all algorithmic components to the framework. The library includes a broad collection of built-in oracles (classical, unified, Pareto-optimal, split cuts, and problem-specific variants) and multiple environment controllers (sequential, in-out stabilized, branch-and-bound). Users can easily integrate custom oracles or environments to explore new algorithmic ideas.

BendersX.jl also provides a growing suite of benchmarking examples—including facility location variants and network interdiction models—designed to support reproducible computational studies and easy comparison across algorithms.

Whether you are developing new Benders decomposition techniques, testing Benders variants, or building scalable optimization applications, BendersX.jl offers a clean, flexible, and extensible platform for working with Benders decomposition in Julia.


## Architecture 
BendersX.jl is organized into three main components, each responsible for part of the algorithm.

### Master
- Represents the master problem in a Benders decomposition.
- Generates candidate points and handles integration of new cuts into the model.
- Provides hooks for user-defined modeling code:
```julia
master = Master(data; customize = customize_master_model!)
```
Master model is given via user-defined function `customize_master_model!`[link]. 

### Oracle
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

## Quick Start
A minimal working example:
```julia
using BendersX
using JuMP

# 1. User-defined data
data = MyData(...)

# 2. Create master model
master = Master(data; customize = customize_master_model!)

# 3. Select oracle and provide subproblem customization
oracle = ClassicalOracle(data, master; customize = customize_sub_model!)

# 4. Choose environment (sequential by default)
env = BendersSeq(master, oracle)

# 5. Solve
log = solve!(env)

```
### Where to go next
- Data structures: How to prepare user-defined input data
- Master customization: How to define variables & constraints
- Oracle customization: How subproblems should be formulated
- Environment behavior: How iterations are controlled

## 4. User Guide
This section explains the key components users must implement when building a Benders model.

### 4.1 Data Structures
Users provide a custom struct to store all instance-specific information:
```julia
struct MyData
    cost::Vector{Float64}
    capacity::Vector{Float64}
    demand::Vector{Float64}
    # …
end
```
This struct is passed into both the master and the oracle, ensuring consistent access to parameters.

### 4.2 Master Model Customization
Users supply a function that defines the JuMP model:
```julia 
function customize_master_model!(model::Model, data::MyData)
    @variable(model, x[1:length(data.cost)], Bin)
    @constraint(model, sum(x) ≤ 10)
    # …
end
```
BendersX.jl wraps the model with additional structures for storing cuts and tracking iteration progress.

### 4.3 Oracle (Subproblem) Customization
The oracle receives master decisions each iteration and solves the subproblem:
```julia
function customize_sub_model!(sub::Model, data::MyData, master_vals)
    # Build subproblem with master_vals fixed
end
```
The oracle must also implement a method that extracts dual information and constructs Benders cuts automatically.

All built-in oracles provide templates to simplify this.