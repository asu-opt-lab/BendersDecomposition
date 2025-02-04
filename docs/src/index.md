

# Introduction

Welcome to BendersDecomposition.jl, a package for solving optimization problems using Benders Decomposition.

BendersDecomposition.jl is built on [JuMP](https://jump.dev), so it supports a number of
open-source and commercial solvers, making it a powerful and flexible tool for
Benders Decomposition.

## Installation

Install `BendersDecomposition.jl` as follows:

```julia
julia> import Pkg

julia> Pkg.add("BendersDecomposition")
```

## License

## Why should you use BendersDecomposition.jl?

Until recently, there were no open-source, generic implementations of the Benders Decomposition algorithm available in the public domain. As a result, practitioners had to develop their own implementations in various languages and styles.

`BendersDecomposition.jl` is designed to be a flexible and extensible package for solving optimization problems using Benders Decomposition. For the beginners, we provide a simple interface to build the model and solve the problem. For the researchers and industry practitioners, we provide the `BendersOracle` interface to customize the modeling, solving and cutting process. In general, the design philosophy of `BendersDecomposition.jl` prioritizes flexibility and efficiency, aiming to deliver a feature-rich and user-friendly experience. 

## Getting started

- Learn the basics of [JuMP](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/getting_started_with_JuMP/) and [Julia](https://jump.dev/JuMP.jl/stable/tutorials/getting_started/getting_started_with_julia/) in the [JuMP documentation](https://jump.dev/JuMP.jl/stable/)
- Follow the tutorials in this manual

If you need help, please open a GitHub issue.

## Citing `BendersDecomposition.jl`

If you use `BendersDecomposition.jl` in your work, please cite the following paper:
