```@meta
CurrentModule = BendersX
```

# Experiments and Reproducibility

The `experiments/` directory is the repository's reproducibility layer. It is
not the primary user tutorial path; instead, it shows how the built-in problem
helpers and algorithmic variants are exercised in benchmark runs.

## Experiment suites

| Suite | Folder or file | Main focus |
| --- | --- | --- |
| Sequential typical | `experiments/1_test_sequential_typical/` | [`BendersSeq`](@ref) with classical, knapsack, and separable oracles |
| Sequential in-out | `experiments/2_test_sequential_in_out_typical/` | [`BendersSeqInOut`](@ref) on the same problem families |
| Sequential disjunctive | `experiments/3_test_sequential_disjunctive/` | [`SplitOracle`](@ref) with sequential execution |
| Callback typical | `experiments/4_test_callback_typical/` | [`BendersBnB`](@ref) with lazy callbacks and root preprocessing variants |
| Callback disjunctive | `experiments/5_test_callback_disjunctive/` | [`BendersBnB`](@ref) with split cuts and user callbacks |
| Specialized sequential | `experiments/6_test_specialized_sequential.jl` | [`SpecializedBendersSeq`](@ref) for split-cut workflows |

Most suite directories contain one script per built-in problem family:

- `ufl.jl`
- `cfl.jl`
- `scfl.jl`
- `snip.jl`

## Reference objectives

`experiments/reference_objectives/` stores benchmark objective values used by
the test suites. The companion script
`experiments/reference_objectives/generate_reference_objectives.jl` regenerates
those reference files.

## Standalone experiment scripts

`experiments/scripts/` contains more targeted experiment entry points, for
example:

- direct MIP baselines;
- CFLP and UFLP benchmarking scripts;
- callback experiments with different disjunctive settings.

Use these files when you want a starting point for your own computational runs
rather than a minimal package tutorial.

## Running the full experiment set

From the repository root:

```bash
julia --project=. experiments/runtests.jl
```

Or from the Julia REPL:

```julia
include("experiments/runtests.jl")
```

## Notes on solver availability

- Sequential experiment suites can often be adapted to open-source solvers.
- Callback-based suites depend on callback-capable solvers.
- Some scripts under `experiments/` are intentionally configured for the
  authors' benchmarking environment and may require local solver or license
  changes before reuse.
