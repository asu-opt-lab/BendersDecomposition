```@meta
CurrentModule = BendersX
```

# [Oracle Guide](@id oracle-guide)

Oracles decide how a candidate master point is evaluated and how new cuts are
generated. In BendersX, every oracle is a subtype of [`AbstractOracle`](@ref)
and implements [`generate_cuts`](@ref).

## Typical LP-based oracles

[`ClassicalOracle`](@ref), [`UnifiedOracle`](@ref), and
[`ParetoOracle`](@ref) all assume that the subproblem built by
[`customize_sub_model!`](@ref) is LP-compatible.

Use them when:

- the recourse model is linear and continuous;
- you want a general-purpose oracle that is independent of problem structure;
- you want to compare classical, unified, and Pareto cuts on the same model.

### Choosing among them

| Oracle | When to prefer it |
| --- | --- |
| [`ClassicalOracle`](@ref) | Baseline implementation and most general starting point |
| [`UnifiedOracle`](@ref) | Unified feasibility/optimality treatment following Fischetti et al. |
| [`ParetoOracle`](@ref) | Stronger cuts when a meaningful core point is available |

### Parameters

- [`ClassicalOracleParam`](@ref) is an alias of [`BasicOracleParam`](@ref).
- [`UnifiedOracleParam`](@ref) adds the `w0` weight for the unified objective
  bound constraint.
- [`ParetoOracleParam`](@ref) requires a problem-specific `core_point`.

```julia
param = ClassicalOracleParam(rtol = 1e-6, atol = 1e-8)
oracle = ClassicalOracle(data, master; customize = customize_sub_model!, param = param)
```

## Problem-specific oracles

BendersX ships problem-specific separation routines for built-in facility
location models.

### [`CFLKnapsackOracle`](@ref)

Use this oracle for capacitated facility location models when the subproblem
matches the CFLP structure. It still builds a model-based subproblem, but it
uses facility-wise knapsack calculations to strengthen the resulting cut.

```julia
oracle = CFLKnapsackOracle(data, master; customize = customize_sub_model!)
```

### [`UFLKnapsackOracle`](@ref)

Use this oracle for uncapacitated facility location. Unlike the model-based
typical oracles, it works directly from [`UFLPData`](@ref) and does not need a
master object to recover a closed-form separation routine.

```julia
param = UFLKnapsackOracleParam(slim = true, add_only_violated_cuts = true)
oracle = UFLKnapsackOracle(data; param = param)
```

## Multi-scenario and separable recourse

[`SeparableOracle`](@ref) wraps a typical oracle prototype and instantiates one
sub-oracle per scenario or block.

It is the standard choice when:

- `t` is vector-valued, one entry per scenario;
- the recourse problem decomposes into independent LP subproblems;
- you want to reuse an existing typical oracle across all scenarios.

```julia
oracle = SeparableOracle(
    data,
    master,
    ClassicalOracle(),
    data.n_scenarios;
    customize = customize_sub_model!,
    sub_oracle_param = ClassicalOracleParam(rtol = 1e-6),
)
```

This pattern is especially useful for [`SCFLPData`](@ref) and
[`SNIPData`](@ref), where scenario-specific subproblems are common.

## Disjunctive separation with [`SplitOracle`](@ref)

[`SplitOracle`](@ref) combines two typical oracles and a
[`DcglpParam`](@ref) configuration to generate split cuts through a disjunctive
cut-generating LP.

Use it when:

- the master problem is mixed-integer;
- you want split cuts in addition to or instead of classical Benders cuts;
- you are experimenting with callback-based or specialized split-cut workflows.

```julia
oracle_kappa = ClassicalOracle(data, master; customize = customize_sub_model!)
oracle_nu = ClassicalOracle(data, master; customize = customize_sub_model!)

dcglp_optimizer = optimizer_with_attributes(CPLEX.Optimizer, MOI.Silent() => true)
dcglp_param = DcglpParam(dcglp_optimizer; verbose = false)
split_param = SplitOracleParam(
    dcglp_param;
    strengthened = true,
    lift = true,
)

oracle = SplitOracle(master, [oracle_kappa, oracle_nu], split_param)
```

Important `SplitOracleParam` controls include:

- `split_index_selection_rule` for choosing the branching disjunction;
- `disjunctive_cut_append_rule` for reusing previously generated split cuts;
- `add_benders_cuts_to_master` for deciding how byproduct Benders cuts are
  passed back to the master;
- `norm`, `reuse_dcglp`, `lift`, and `strengthened` for DCGLP behavior.

## Practical guidance

- Start with [`ClassicalOracle`](@ref) unless you have a clear reason to use a
  stronger or problem-specific oracle.
- Prefer [`SeparableOracle`](@ref) whenever the recourse structure is naturally
  scenario-wise and `t` is vector-valued.
- Use [`CFLKnapsackOracle`](@ref) and [`UFLKnapsackOracle`](@ref) only when the
  built-in facility-location structure matches your model.
- Use [`SplitOracle`](@ref) together with
  [Environment Guide](@ref environment-guide) when you need callback-based or
  specialized split-cut workflows.
