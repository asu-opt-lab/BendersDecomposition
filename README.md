# BendersX.jl

[![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://asu-opt-lab.github.io/BendersX.jl/)
[![Julia](https://img.shields.io/badge/julia-v1.10.4-blue.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

*A modular, hierarchical, and plug-and-play framework for Benders decomposition in Julia.*
BendersX.jl provides a fully extensible architecture for implementing, experimenting with, and benchmarking Benders decomposition algorithms. Unlike problem-specific or monolithic implementations, BendersX.jl separates the algorithm into independent components—**Master**, **Oracle**, and **Environment**—each of which can be independently customized or replaced. This architecture enables rapid prototyping of new algorithmic ideas, fair comparison across Benders variants, and reproducible computational studies.

## Key Features
- **Extensive Oracle Library**
    Includes built-in oracles that generate:
    - Classical Benders cuts (optimality & feasibility cuts)
    - Unified Benders cuts
    - Pareto-optimal cuts
    - Problem-specific Benders cuts (e.g., knapsack-based for facility location problems)
    - Split cuts for Benders reformulation
- **Flexible Environment Library**
    Support multiple environments that implement:
    - Sequential Benders
    - In-out stabilized Benders
    - Branch-and-bound Benders
- **JuMP-Based Modeling**
    Users write master and subproblem models using standard JuMP syntax; BendersX.jl handles all algorithmic logics behind the scenes.
- **Plug-and-Play Extensibility**
    Oracles and environments can be swapped freely with no changes to modeling code. Users can also implement custom oracles or environments and integrate them directly into the framework.
- **Benchmarking Suite**
    Includes ready-to-run examples for:
    - (Uncapacitated, Capacitated, Stochastic Capacitated) Facility Location Problem (UFLP, CFLP, SCFLP)
    - Stochastic Network Interdiction Problem (SNIP)

## Installation
To set up the project:
```julia
"""
Installation instructions will be added after package registration.
"""
```

## Quick Start
```julia
using BendersX
using JuMP

# 1. include here user-defined data struct
data = MyData(...)
# 2. include here master and subproblem customization
master = Master(data; customize = customize_master_model!)
oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
env = BendersSeq(master, oracle)
log = solve!(env)
```
1. For details on **data structures**, see:
➡️ [Documentation: Data Guide](#) <!-- replace # with actual URL -->
2. For details on **model customization functions**, see:
➡️ [Documentation: Model Customization Guide](#) <!-- replace # with actual URL -->

For the **full documentation homepage**, visit:
➡️ [Documentation](#)

## Examples
Several Julia scripts are provided to run Benders algorithms on different problem instances.
See the `experiments/` directory for more information.

## Contributing
This repository is under active development, and contributions are welcome.
Please submit issues for bugs or feature requests, and open a pull request for proposed changes.
For major updates, we recommend discussing your idea in an issue before submitting code.

## License
Copyright © 2025 Arizona State University.
Released under the MIT License (see [LICENSE](LICENSE) file for details).







