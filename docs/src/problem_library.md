```@meta
CurrentModule = BendersX
```

# [Problem Library](@id problem-library)

BendersX includes several built-in data containers, readers, and modeling
helpers for benchmark problem families. These helpers are intended both as
ready-to-use examples and as templates for your own models.

## Artifact-backed datasets

The built-in readers use Julia artifacts. On first use, the corresponding
dataset is downloaded from the repository's release assets and cached locally.

- The default dataset location is resolved internally from the package
  artifacts.
- Every reader also accepts an override keyword such as `filepath` or
  `base_dir` so you can point it at a local directory instead of the packaged
  artifact.

## Built-in problem families

| Problem family | Data type | Main readers | Typical oracle patterns |
| --- | --- | --- | --- |
| Capacitated facility location | [`CFLPData`](@ref) | [`read_GK_data`](@ref), [`read_cflp_benchmark_data`](@ref), [`read_cfl_file`](@ref) | [`ClassicalOracle`](@ref), [`UnifiedOracle`](@ref), [`ParetoOracle`](@ref), [`CFLKnapsackOracle`](@ref) |
| Uncapacitated facility location | [`UFLPData`](@ref) | [`read_uflp_benchmark_data`](@ref), [`read_Simple_data`](@ref) | [`ClassicalOracle`](@ref), [`UFLKnapsackOracle`](@ref) |
| Stochastic capacitated facility location | [`SCFLPData`](@ref) | [`read_stochastic_capacited_facility_location_problem`](@ref) | [`SeparableOracle`](@ref) over a typical oracle |
| Stochastic network interdiction | [`SNIPData`](@ref) | [`read_snip_data`](@ref) | [`SeparableOracle`](@ref) or custom research workflows |

## CFLP

[`CFLPData`](@ref) stores facility capacities, customer demands, fixed costs,
and the transportation cost matrix.

Available readers:

- [`read_GK_data`](@ref) for JSON-form random instances;
- [`read_cflp_benchmark_data`](@ref) for benchmark files from the packaged
  `cflp_locssall` artifact;
- [`read_cfl_file`](@ref) for `.cfl` formatted instances.

Typical workflows:

- benchmark-friendly sequential runs with [`ClassicalOracle`](@ref);
- specialized separation with [`CFLKnapsackOracle`](@ref);
- direct monolithic baselines through [`customize_mip_model!`](@ref).

## UFLP

[`UFLPData`](@ref) stores facility-opening costs, customer demands, and the
cost matrix for uncapacitated facility location.

Readers:

- [`read_uflp_benchmark_data`](@ref) for packaged benchmark files;
- [`read_Simple_data`](@ref) for the packaged "Simple" instances.

The closed-form [`UFLKnapsackOracle`](@ref) is the main problem-specific oracle
for this family.

## SCFLP

[`SCFLPData`](@ref) stores a vector of demand scenarios and therefore usually
pairs with vector-valued `t` variables in the master problem.

The most common workflow is:

1. build the master with one `t` entry per scenario;
2. implement [`customize_sub_model!`](@ref) with a `scen_idx` argument;
3. wrap a typical oracle prototype in [`SeparableOracle`](@ref).

## SNIP

[`SNIPData`](@ref) represents stochastic network interdiction instances with
scenario probabilities, arc data, and a budget. As with SCFLP, it naturally
fits a separable multi-scenario workflow.

[`read_snip_data`](@ref) expects:

- `instance_no`;
- `snip_no`;
- `budget`.

## Built-in MIP helpers

For the shipped problem data types, the repository provides
[`customize_mip_model!`](@ref) methods that build direct MIP baselines. These
methods are useful for experiment scripts and regression checks, especially
when comparing Benders against a monolithic solve.

## What is intentionally not documented here

The repository contains some experimental or partially wired assets outside the
public package surface. This page only documents problem families that are
currently exported or marked `public` from `BendersX`.
