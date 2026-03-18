```@meta
CurrentModule = BendersX
```

# BendersX.jl

**BendersX.jl** is a Julia framework for building, comparing, and extending
Benders decomposition algorithms with interchangeable **masters**,
**oracles**, and **execution environments**.

The package is designed for two common workflows:

- building a new decomposition model from JuMP components;
- running built-in facility-location or network-interdiction benchmarks with
  different oracle and environment combinations.

## Core workflow

Every BendersX solve follows the same assembly pattern:

```julia
using BendersX

data = MyData(...)
master = Master(data; customize = customize_master_model!)
oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
env = BendersSeq(master, oracle)
log = solve!(env)
```

The modeling details come from user or problem-specific implementations of
[`customize_master_model!`](@ref) and [`customize_sub_model!`](@ref). The
algorithmic behavior comes from the selected oracle and environment.

## Installation

Install from GitHub:

```julia
import Pkg
Pkg.add(url = "https://github.com/asu-opt-lab/BendersX.jl")
```

For local development from a clone:

```julia
import Pkg
Pkg.develop(path = "/path/to/BendersDecomposition")
```

## Solver support

- Sequential workflows such as [`BendersSeq`](@ref) and
  [`BendersSeqInOut`](@ref) can be demonstrated with open-source LP/MIP
  solvers such as `HiGHS`.
- Callback-based workflows such as [`BendersBnB`](@ref),
  [`LazyCallback`](@ref), and [`UserCallback`](@ref) require a solver with
  callback support.
- The built-in problem helper methods in `src/problems/` currently configure
  commercial solvers in source code. Treat them as reference implementations
  and adapt the optimizer choice if your local setup differs.

## What is documented here

- [Getting Started](@ref getting-started) gives the shortest sequential example that runs with
  a small built-in data object and `HiGHS`.
- [Architecture](@ref architecture) explains how masters, oracles, and
  environments fit together.
- [Modeling Guide](@ref modeling-guide) documents the user-facing model
  customization contracts, including generalized bound constraints.
- [Problem Library](@ref problem-library) covers the built-in data types,
  readers, and artifact-backed benchmark sets.
- [Oracle Guide](@ref oracle-guide) and [Environment Guide](@ref environment-guide)
  explain when to use each algorithmic component.
- [Extending BendersX](@ref extending-bendersx) covers custom oracle and
  environment implementations.
- [API Reference](@ref api) groups the public surface by workflow area.

## Built-in algorithm families

| Family | Primary entry points | Typical use |
| --- | --- | --- |
| Sequential | [`BendersSeq`](@ref), [`BendersSeqInOut`](@ref) | Standard cutting-plane workflows and debugging |
| Callback-based | [`BendersBnB`](@ref), [`LazyCallback`](@ref), [`UserCallback`](@ref) | Mixed-integer masters with solver callbacks |
| Problem-specific oracles | [`CFLKnapsackOracle`](@ref), [`UFLKnapsackOracle`](@ref) | Specialized facility-location separation |
| Disjunctive workflows | [`SplitOracle`](@ref), [`SpecializedBendersSeq`](@ref) | Split cuts and research variants |

!!! info "Complete example"
    For a benchmark-backed example built from the shipped CFLP helpers, see
    [CFLP Demo](@ref cflp-demo).
