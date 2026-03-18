```@meta
CurrentModule = BendersX
```

# [Environment Guide](@id environment-guide)

Environments determine **how** a master and oracle interact. They own the solve
loop, stopping rules, callback registration, and logging behavior.

## Sequential environments

### [`BendersSeq`](@ref)

[`BendersSeq`](@ref) is the default cutting-plane environment:

- one master solve per iteration;
- one oracle call per iteration;
- best for initial model validation and baseline comparisons.

```julia
seq_param = BendersSeqParam(time_limit = 600.0, gap_tolerance = 1e-4, verbose = true)
env = BendersSeq(master, oracle; param = seq_param)
log = solve!(env)
```

### [`BendersSeqInOut`](@ref)

[`BendersSeqInOut`](@ref) adds in-out stabilization. Use it when the plain
sequential loop oscillates or progresses slowly and you can provide a sensible
stabilizing point.

```julia
inout_param = BendersSeqInOutParam(
    time_limit = 600.0,
    gap_tolerance = 1e-4,
    stabilizing_x = fill(0.5, master.dim_x),
    α = 0.9,
    λ = 0.1,
)
env = BendersSeqInOut(master, oracle; param = inout_param)
```

## Callback-based branch and bound

[`BendersBnB`](@ref) integrates Benders cuts into a MIP solver through lazy and
user callbacks.

Use it when:

- the master is mixed-integer;
- your solver supports callbacks;
- you want to add cuts at integer or fractional nodes instead of solving a
  standalone sequential loop.

```julia
bnb_param = BendersBnBParam(time_limit = 1800.0, gap_tolerance = 1e-4, verbose = true)
env = BendersBnB(master, oracle; param = bnb_param)
log = solve!(env)
```

### Callback components

- [`LazyCallback`](@ref) adds cuts at integer nodes and is the default callback
  wrapper for typical oracles.
- [`UserCallback`](@ref) adds cuts at fractional nodes, usually for stronger
  disjunctive separation.
- [`NoUserCallback`](@ref) disables user cuts.

## Root preprocessing

[`BendersBnB`](@ref) supports optional root-node preprocessing before callback
search begins.

### No preprocessing

```julia
root_preprocessing = NoRootNodePreprocessing()
env = BendersBnB(master, root_preprocessing, LazyCallback(oracle), NoUserCallback(); param = bnb_param)
```

### Sequential preprocessing

[`RootNodePreprocessing`](@ref) runs a sequential environment on the LP-relaxed
master to add initial cuts:

```julia
root_preprocessing = RootNodePreprocessing(
    oracle,
    BendersSeq,
    BendersSeqParam(time_limit = 120.0, gap_tolerance = 1e-6, verbose = false),
)
env = BendersBnB(master, root_preprocessing, LazyCallback(oracle), NoUserCallback(); param = bnb_param)
```

You can swap `BendersSeq` for [`BendersSeqInOut`](@ref) when stabilization at
the root node is useful.

### Disjunctive root preprocessing

[`DisjunctiveRootNodePreprocessing`](@ref) runs a typical phase followed by a
disjunctive phase. It is intended for split-cut workflows.

## Combining lazy and user callbacks

When you want both standard Benders cuts at integer nodes and stronger split
cuts at fractional nodes, use a typical oracle in the lazy callback and a
[`SplitOracle`](@ref) inside [`UserCallback`](@ref).

```julia
kappa = ClassicalOracle(data, master; customize = customize_sub_model!)
nu = ClassicalOracle(data, master; customize = customize_sub_model!)
dcglp_param = DcglpParam(optimizer_with_attributes(CPLEX.Optimizer, MOI.Silent() => true); verbose = false)
split_param = SplitOracleParam(dcglp_param; strengthened = true, lift = true)

lazy_callback = LazyCallback(kappa)
user_callback = UserCallback(SplitOracle(master, [kappa, nu], split_param); params = UserCallbackParam(frequency = 1))

env = BendersBnB(master, NoRootNodePreprocessing(), lazy_callback, user_callback; param = bnb_param)
```

## Specialized split-cut workflow

[`SpecializedBendersSeq`](@ref) is a research-oriented sequential environment
for split-cut workflows. It expects a [`SplitOracle`](@ref) configured with the
rules required by the implementation:

- `split_index_selection_rule = LargestFractional()`
- `disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices()`

```julia
split_param = SplitOracleParam(
    dcglp_param;
    split_index_selection_rule = LargestFractional(),
    disjunctive_cut_append_rule = DisjunctiveCutsSmallerIndices(),
)
split_oracle = SplitOracle(master, [kappa, nu], split_param)
env = SpecializedBendersSeq(master, split_oracle)
```

## Choosing an environment

| Environment | Best first use |
| --- | --- |
| [`BendersSeq`](@ref) | Standard baseline and debugging |
| [`BendersSeqInOut`](@ref) | Sequential runs that benefit from stabilization |
| [`BendersBnB`](@ref) | Mixed-integer masters with callback-capable solvers |
| [`SpecializedBendersSeq`](@ref) | Split-cut experiments and specialized research workflows |
