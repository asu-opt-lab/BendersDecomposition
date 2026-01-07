# BendersX.jl

Welcome to **BendersX.jl** — a modular, plug-and-play framework for Benders decomposition algorithms in Julia.

[![License](https://img.shields.io/badge/license-MIT-lightgrey.svg)](LICENSE)
[![Docs](https://img.shields.io/badge/docs-source-blue.svg)](docs/)


---

## Introduction

**BendersX.jl** is a modular and extensible framework for Benders decomposition in Julia. It supports both standard implementations and experimental extensions within a unified, principled design.

The **“X”** in BendersX signals extension and exploration: new application domains and methodological extensions beyond classical Benders decomposition — including alternative cut-generation strategies, stabilization mechanisms, and execution logic. The framework’s plug-and-play architecture makes it easy to implement, combine, and evaluate these techniques with minimal boilerplate.

---

> **Quick overview**
>
> - **Components:** Master · Oracle (cut generator) · Environment (execution controller)  
> - **Architecture:** JuMP-based modeling with a modular, hierarchical algorithmic framework  
> - **Design goal:** Plug-and-play extensibility for rapid prototyping and reproducible experimentation

---

## Quick start

### Installation

Install BendersX with the Julia package manager. At the `julia>` prompt:

```julia
import Pkg
Pkg.add("BendersX")
```

If you are developing the package locally, activate the package project and `dev` it:

> ```julia
> julia> ]
> pkg> dev .
> ```

### Minimal workflow

This example shows a minimal end-to-end workflow. 

```julia
using BendersX, JuMP

# 1. User-defined data
data = MyData(...)

# 2. Create master and provide master problem customization
master = Master(data; customize = customize_master_model!)

# 3. Select oracle and provide subproblem customization
oracle = ClassicalOracle(data, master; customize = customize_sub_model!)

# 4. Choose environment
env = BendersSeq(master, oracle)

# 5. Solve
log = solve!(env)
```

#### Modeling

Users provide master and subproblem formulations through *customization functions* written in standard JuMP syntax. See the package docs for concrete examples and the full modeling interface. If you are new to JuMP, consult the JuMP documentation: [https://jump.dev/JuMP.jl/stable/](https://jump.dev/JuMP.jl/stable/).

---

#### Built-in variants

#### Oracle variants (examples)

| Oracle type         | Description                                                          |
| ------------------- | -------------------------------------------------------------------- |
| `ClassicalOracle`   | Standard Benders cut generation based on dual information            |
| `UnifiedOracle`     | Unified handling of feasibility and optimality cuts                  |
| `ParetoOracle`      | Produces Pareto-optimal Benders cuts                                 |
| `SeparableOracle`   | Wrapper oracle for separable subproblems      |
| `UFLKnapsackOracle` | Knapsack-based oracle for uncapacitated facility location            |
| `CFLKnapsackOracle` | Knapsack-based oracle for capacitated facility location              |
| `SplitOracle`       | Generates split cuts to strengthen the master relaxation |

#### Environment variants (examples)

| Environment type  | Execution strategy                            |
| ----------------- | --------------------------------------------- |
| `BendersSeq`      | Classical sequential Benders decomposition    |
| `BendersSeqInOut` | Sequential Benders with in–out stabilization  |
| `BendersBnB`      | Branch-and-bound integrated with Benders cuts |

> These are representative built-ins. See the documentation for the full list and configuration options.

---

## Key features

* **Plug-and-play architecture:** swap or extend Master, Oracle, and Environment components without rewriting models.
* **JuMP-native modeling:** define master and subproblems using standard JuMP expressions and containers.
* **Algorithmic variants:** support for classical, unified, Pareto, split cuts, in–out stabilization, B&B integration, and more.
* **Benchmark-ready:** included examples and suites for facility location, network interdiction, and other problems to support reproducible experiments.

---

## Next steps

- Refer to the **Tutorials** documentation for worked examples and step-by-step guides.
- Refer to the **User Guide** for advanced usage patterns and customization guidance.
- Consult the **API** documentation for detailed descriptions of the Master, Oracle, and Environment interfaces.

---

## Contributing

Contributions, bug reports, and enhancements are welcome. Please open issues or pull requests on GitHub. Follow the repository's coding guidelines and run the test suite locally before submitting PRs.

---

## License

BendersX.jl is released under the **MIT License**. See [LICENSE](https://github.com/asu-opt-lab/BendersX.jl/blob/main/LICENSE) for details.
