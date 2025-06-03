# Disjunctive Benders Decomposition

[![Julia](https://img.shields.io/badge/julia-v1.10.4-blue.svg)](https://julialang.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

A Julia implementation of disjunctive Benders decomposition algorithms for solving mixed-integer programming problems, developed as part of research on "Disjunctive Benders Decomposition".

## Overview

This repository contains the source code and computational experiments for disjunctive Benders decomposition methods. The implementation extends classical Benders decomposition by incorporating disjunctive cuts to improve convergence and solution quality for mixed-integer programming problems.

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
- `Dcglp`: Disjunctive Cut Generating Linear Program
- `SpecializedBendersSeq`: Specialized sequential implementation

### Oracle Types
- `ClassicalOracle`: Traditional Benders subproblem oracle
- `KnapsackOracle`: Knapsack technique based oracle
- `DisjunctiveOracle`: Disjunctive programming-based oracle
- `SeparableOracle`: Oracle for separable subproblems

## Problem Examples

The `example/` directory contains implementations for several classic optimization problems:

- **UFLP**: Uncapacitated Facility Location Problem
- **CFLP**: Capacitated Facility Location Problem  
- **SCFLP**: Stochastic Facility Location Problem
- **MCNDP**: Multi-Commodity Network Design Problem
- **SNIP**: Stochastic Network Interdiction Problem

## Installation

To set up the project:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

## Usage

We provide several scripts to run the algorithms on different problem instances. Please refer to the `scripts/` directory for more details.

## Testing

Run the test suite:
```bash
julia test/runtests.jl
```

Or run specific tests:
```bash
./test/runtests.sh
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
│   ├── algorithms/          # Core decomposition algorithms
│   ├── modules/            # Oracle implementations and components
│   ├── utils/              # Utility functions and helpers
│   └── types.jl            # Type definitions and exports
├── test/                   # Comprehensive test suite
├── example/                # Problem-specific implementations
│   ├── uflp/              # Uncapacitated facility location
│   ├── cflp/              # Capacitated facility location
│   ├── scflp/             # Single-commodity flow location
│   ├── mcndp/             # Multi-commodity network design
│   └── snip/              # Stochastic network interdiction
└── Project.toml           # Julia project configuration
```

## Contributing

We welcome contributions! Please feel free to submit issues and pull requests. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.








