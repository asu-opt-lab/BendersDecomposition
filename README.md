# BendersX.jl

Welcome to **BendersX.jl** — a modular, plug-and-play framework for Benders decomposition algorithms in Julia.

[![Julia](https://img.shields.io/badge/julia-v1.11%2B-blue.svg)](https://julialang.org/)
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
> - **Repository layout:** package source in `src/`, docs in `docs/`, and benchmark scripts in `experiments/`

---

## Installation

Install BendersX with the Julia package manager.

To add the package from GitHub:
```julia
import Pkg
Pkg.add(url = "https://github.com/asu-opt-lab/BendersDecomposition.git")
```

If you are developing the repository locally:
```julia
import Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Quick start

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

The default `using BendersX` workflow brings in the core constructors and `solve!`. Advanced extension hooks and problem-specific helpers remain public APIs and can be accessed with either `import BendersX: ...` or explicit qualification such as `BendersX.read_cflp_benchmark_data`.

#### Built-in variants

#### Oracle variants (examples)

| Oracle type         | Description                                                          |
| ------------------- | -------------------------------------------------------------------- |
| `ClassicalOracle`   | Standard Benders cut generation based on dual information            |
| `UnifiedOracle`     | Unified handling of feasibility and optimality cuts                  |
| `ParetoOracle`      | Produces Pareto-optimal Benders cuts                                 |
| `SeparableOracle`   | Wrapper oracle for separable subproblems                             |
| `UFLKnapsackOracle` | Knapsack-based oracle for uncapacitated facility location            |
| `CFLKnapsackOracle` | Knapsack-based oracle for capacitated facility location              |
| `SplitOracle`       | Generates split cuts to strengthen the master relaxation             |

#### Environment variants (examples)

| Environment type        | Execution strategy                                |
| ----------------------- | ------------------------------------------------- |
| `BendersSeq`            | Classical sequential Benders decomposition        |
| `BendersSeqInOut`       | Sequential Benders with in-out stabilization      |
| `SpecializedBendersSeq` | Specialized sequential workflow for split cuts    |
| `BendersBnB`            | Branch-and-bound integrated with Benders cuts     |

> These are representative built-ins. See the documentation for the full list and configuration options.

---

## Repository layout

This repository contains the package source code, benchmark loaders, documentation, and computational experiments for `BendersX`.

- `src/BendersX.jl` centralizes the package public API and exports
- `src/modules/` contains masters, environments, oracles, and callback logic
- `src/problems/` contains benchmark readers, problem models, and specialized problem-specific oracles
- `src/artifact_utils.jl` and `Artifacts.toml` manage artifact-backed benchmark datasets downloaded on demand
- `docs/` contains the standalone Documenter site, tutorials, and API pages
- `experiments/` contains experiment scripts and reference objectives used for benchmarking
- `test/` contains public API and unit tests

---

## Key features

- **Plug-and-play architecture:** swap or extend Master, Oracle, and Environment components without rewriting models.
- **JuMP-native modeling:** define master and subproblems using standard JuMP expressions and containers.
- **Algorithmic variants:** support for classical, unified, Pareto, split cuts, in-out stabilization, branch-and-bound integration, and more.
- **Artifact-backed benchmark data:** built-in readers for CFLP, UFLP, SCFLP, and SNIP obtain packaged datasets lazily when needed.
- **Benchmark-ready workflow:** tutorials, docs, and experiment suites are organized for reproducible experimentation.

---

## Testing and docs

Run the test suite:

```bash
julia --project=. test/runtests.jl
```

Build the documentation locally:

```bash
julia --project=docs docs/make.jl
```

Benchmark and algorithm comparison scripts live under `experiments/`, including sequential, in-out, callback, and disjunctive experiment setups.

---

## Next steps

- Refer to the **Tutorials** documentation for worked examples and step-by-step guides.
- Refer to the **User Guide** for advanced usage patterns and customization guidance.
- Consult the **API** documentation for detailed descriptions of the Master, Oracle, and Environment interfaces.

## Contributing

Contributions, bug reports, and enhancements are welcome. Please open issues or pull requests on GitHub. Follow the repository's coding guidelines and run the test suite locally before submitting PRs.

---

## License

BendersX.jl is released under the **MIT License**. See [LICENSE](LICENSE) for details.
