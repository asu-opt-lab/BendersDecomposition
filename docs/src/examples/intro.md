```@meta
EditURL = "intro.jl"
```

# Introduction

In this package, we provide some instances and some cut strategies.

```mermaid
graph TD
  A[Instances] --> B[Capacitated Facility Location Problem]
  A --> C[Uncapacitated Facility Location Problem]
  A --> D[Stochastic Capacitated Facility Location Problem]
  A --> E[Multicommodity Capacitated Network Design Problem]
  A --> F[Stochastic Network Interdiction Problem]
  B --> G[Classical Cut]
  C --> G
  D --> G
  E --> G
  F --> G
  B --> H[Knapsack Cut]
  C --> H
  D --> H
  E --> H
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

