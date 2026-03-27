# BendersX.jl

[![Julia](https://img.shields.io/badge/julia-v1.11%2B-blue.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

A Julia framework for Benders decomposition research and experimentation, with modular Master, Oracle, and Environment components and built-in support for classical, unified, Pareto, separable, and disjunctive variants.

## Overview

This repository contains the source code, benchmark problem loaders, tutorials, and computational experiments for the `BendersX` package. The codebase is organized as a single Julia package rooted at `src/BendersX.jl`, with package internals grouped under `src/modules/`, `src/problems/`, and `src/utils/`.

## Key Features

- **Multiple Algorithm Variants**: Implementation of sequential and callback-based Benders decomposition algorithms
- **Disjunctive Cuts**: Integration of disjunctive programming techniques for enhanced cut generation
- **Flexible Oracle System**: Modular oracle design supporting different subproblem types (typical, disjunctive, separable)
- **Comprehensive Testing**: Extensive test suite with multiple problem instances
- **Multiple Problem Types**: Support for facility location problems (UFLP, CFLP, SCFLP), network design (MCNDP), and other optimization problems

## Algorithm Implementations

### Core Algorithms
- `BendersSeq`: Sequential Benders decomposition
- `BendersSeqInOut`: Sequential variant with in-out technique
- `BendersBnB`: Branch-and-bound Benders decomposition  
- `Dcglp`: Split Cut Generating Linear Program
- `SpecializedBendersSeq`: Specialized sequential implementation

### Oracle Types
- `ClassicalOracle`: Traditional Benders subproblem oracle
- `UnifiedOracle`: Unified treatment of feasibility and optimality cuts
- `ParetoOracle`: Pareto-optimal cut generation
- `KnapsackOracle`: Knapsack technique based oracle
- `SplitOracle`: Disjunctive programming-based oracle
- `SeparableOracle`: Oracle for separable subproblems

## Built-in Problem Loaders

The `src/problems/` directory contains benchmark readers and helpers for:
- **UFLP**: Uncapacitated Facility Location Problem
- **CFLP**: Capacitated Facility Location Problem  
- **SCFLP**: Stochastic Capacitated Facility Location Problem
- **MCNDP**: Multi-Commodity Network Design Problem
- **SNIP**: Stochastic Network Interdiction Problem

The `experiments/` directory contains scripts used to benchmark these formulations under different oracle and environment choices.

## Installation

To work from a local clone:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

To add the package from GitHub:
```julia
using Pkg
Pkg.add(url = "https://github.com/asu-opt-lab/BendersX.jl")
```

## Usage

The package entry point is `using BendersX` for the core workflow (`Master`, the standard oracle/environment constructors, and `solve!`). Advanced extension hooks and problem-specific helpers remain public APIs, but may require either `import BendersX: ...` or explicit qualification such as `BendersX.read_cflp_benchmark_data`.

Tutorials and worked examples live under `docs/src/tutorials/`, and benchmark scripts live under `experiments/`.

## Testing

Run the test suite:
```bash
julia --project=. test/runtests.jl
```

Build the documentation locally:
```bash
julia --project=docs docs/make.jl
```

The test suite includes:
1. Sequential typical Benders decomposition
2. Sequential in-out typical Benders decomposition  
3. Sequential disjunctive Benders decomposition
4. Callback typical Benders decomposition
5. Callback disjunctive Benders decomposition
6. Specialized sequential Benders decomposition

## Project Structure

```
├── src/
│   ├── BendersX.jl         # Package entry point and exports
│   ├── modules/            # Masters, environments, oracles, and callbacks
│   ├── problems/           # Benchmark data readers and problem-specific helpers
│   ├── utils/              # Modeling and algorithm utilities
│   ├── artifact_utils.jl   # Artifact-backed data access helpers
│   └── types.jl            # Core type definitions and parameters
├── docs/                   # Standalone Documenter site and tutorials
├── experiments/            # Experiment and benchmarking scripts
├── test/                   # Comprehensive test suite
└── Project.toml           # Julia project configuration
```

## Contributing

This repository is actively under development and we welcome contributions! Feel free to submit issues for bugs or feature requests, and pull requests for code changes. For major modifications, please open an issue first to discuss your proposal. We appreciate all contributions, from bug fixes to documentation improvements.

## License

Copyright © 2025 Arizona State University.
Released under the MIT License (see [LICENSE](LICENSE) file for details).




