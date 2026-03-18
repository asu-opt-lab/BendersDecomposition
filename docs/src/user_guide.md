```@meta
CurrentModule = BendersX
```

# [Architecture](@id architecture)

BendersX is organized around three replaceable components:

- a **master** that stores the JuMP model seen by the algorithm;
- an **oracle** that evaluates a candidate point and produces cuts;
- an **environment** that controls the solve loop or callback workflow.

This separation keeps modeling decisions, cut-generation logic, and execution
strategy independent.

## Component boundaries

### Master

[`Master`](@ref) stores the JuMP model, the flattened coupling variables `x`,
the auxiliary approximation variables `t`, and the objective coefficients used
to evaluate bounds.

Use the master layer when you need to change:

- the first-stage formulation;
- the shape of the coupling variables passed into the subproblem;
- the meaning and dimension of the approximation variables `t`.

### Oracle

An [`AbstractOracle`](@ref) receives a candidate `(x, t)` point and returns a
tuple of the form `(is_in_L, hyperplanes, f_x)` through
[`generate_cuts`](@ref).

Use the oracle layer when you need to change:

- how subproblems are solved;
- what cuts are generated;
- how feasibility and optimality are detected;
- whether the subproblem structure is classical, separable, problem-specific,
  or disjunctive.

### Environment

An [`AbstractBendersEnv`](@ref) owns the algorithmic control flow. Sequential
environments repeatedly solve the master and query an oracle. Callback-based
environments integrate cut generation into a MIP solver. Specialized
environments may impose additional rules on the master and oracle pair.

Use the environment layer when you need to change:

- the iteration logic;
- stabilization and stopping rules;
- callback behavior;
- root preprocessing;
- how logs and progress data are collected.

## Type hierarchies

The package exposes the main hierarchy explicitly so that new implementations
can be built incrementally rather than from scratch.

### Environment hierarchy

![Environment hierarchy](EnvHierarchy.png)

### Oracle hierarchy

![Oracle hierarchy](OracleHierarchy.png)

## How the pieces fit together

1. Build a data container that is a subtype of [`AbstractData`](@ref).
2. Construct a [`Master`](@ref) with [`customize_master_model!`](@ref).
3. Choose an oracle such as [`ClassicalOracle`](@ref),
   [`SeparableOracle`](@ref), or [`SplitOracle`](@ref).
4. Choose an environment such as [`BendersSeq`](@ref),
   [`BendersSeqInOut`](@ref), [`BendersBnB`](@ref), or
   [`SpecializedBendersSeq`](@ref).
5. Call [`solve!`](@ref).

## Choosing an extension point

- If you are changing the first-stage formulation or the coupling variables,
  update the master customization function.
- If you are changing how cuts are derived, add or modify an oracle.
- If you are changing how and when the master and oracle interact, add or
  modify an environment.

The dedicated extension contracts are documented in
[Extending BendersX](@ref extending-bendersx).
