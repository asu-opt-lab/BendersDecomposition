## Summary of Experiment Test Suites

The `experiments/` directory contains a structured collection of runnable test
suites. Each suite focuses on a specific **execution environment**, and within
each suite, experiments are run on multiple benchmark problems:

* Uncapacitated Facility Location Problem (UFLP)
* Capacitated Facility Location Problem (CFLP)
* Stochastic Capacitated Facility Location Problem (SCFLP)
* Stochastic Network Interdiction Problem (SNIP)

| Test suite | Folder | Oracles demonstrated | Environments demonstrated |
| --- | --- | --- | --- |
| **Sequential (typical)** | `1_test_sequential_typical/` | Typical oracles, including `ClassicalOracle`, problem-specific knapsack oracles, and `SeparableOracle` (for SCFLP and SNIP) | `BendersSeq` |
| **Sequential In-Out (typical)** | `2_test_sequential_in_out_typical/` | Typical oracles (same as above) | `BendersSeqInOut` |
| **Sequential (disjunctive)** | `3_test_sequential_disjunctive/` | `DistanceNormOracle` combined with typical oracles | `BendersSeq` |
| **Callback (typical)** | `4_test_callback_typical/` | Typical oracles (classical, knapsack-based, and separable) | `BendersBnB` with lazy callbacks and varying root preprocessing via subtypes of `AbstractBendersSeq` |
| **Callback (disjunctive)** | `5_test_callback_disjunctive/` | `DistanceNormOracle` combined with typical oracles | Same as above, additionally including user callbacks |

---


## Structure of Each Test Suite

Within each test-suite directory (e.g. `1_test_sequential_typical/`), the
corresponding `runtest.jl` typically includes:

```text
ufl.jl    # Uncapacitated Facility Location
cfl.jl    # Capacitated Facility Location
scfl.jl   # Stochastic Capacitated Facility Location
snip.jl   # Stochastic Network Interdiction
```

Each of these files:

* loads problem-specific data, models, and oracles from `src/problems/...`,
* constructs appropriate masters, oracles, and environments,
* runs the algorithm and validates correctness (often against a reference MIP).

---

## How to Run All Experiments

To execute the full experimental test suite:

```bash
julia --project=. experiments/runtests.jl
```

Or from the Julia REPL:

```julia
julia> include("experiments/runtests.jl")
```

You can also run an individual suite or problem by including the corresponding
`runtest.jl` or problem file directly.

---
