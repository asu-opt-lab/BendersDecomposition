```@meta
CurrentModule = BendersX
```

In BendersX, **environments** control *how* the Benders algorithm is executed.
They orchestrate the interaction between the master object and the oracle(s),
determine iteration order, and enforce algorithmic constraints.

This design allows users to swap environments and adjust their behavior without
modifying either the master or oracle objects.

This section explains:

* How to switch between different Benders environments
* How environment parameters and subcomponents affect algorithm execution

---

## Swapping Execution Environments

All execution environments in BendersX are subtypes of
[`AbstractBendersEnv`](@ref). The same master and oracle objects can therefore be
reused across different environments.

For example, switching from a sequential environment to a branch-and-bound
environment only requires changing the environment constructor:

```julia
# Standard sequential Benders
env = BendersSeq(master, oracle)

# Branch-and-bound Benders
env = BendersBnB(master, oracle)
```

The master and oracle objects remain unchanged.

---

## Sequential vs. Branch-and-Bound Environments

BendersX provides two main classes of execution environments.

### Sequential Environments

Sequential environments are subtypes of [`AbstractBendersSeq`](@ref) and execute
the Benders algorithm in a serial, iteration-by-iteration fashion.

Typical use cases include:

* LP master problems
* Continuous relaxations of MILP masters
* Algorithm development and debugging

Example:

```julia
env = BendersSeq(master, oracle)
```

### Branch-and-Bound Environments

Branch-and-bound environments are subtypes of
[`AbstractBendersBnB`](@ref) and integrate Benders decomposition into a
branch-and-bound framework.

These environments are typically used for:

* Integer or mixed-integer master problems
* Disjunctive or logic-based Benders methods
* Callback-based cut generation

Example:

```julia
env = BendersBnB(master, oracle)
```

---

## Adjusting Environment Behavior via Parameters and Subcomponents

Each environment owns a parameter object that controls stopping criteria,
logging, and algorithmic options. In addition, some environments expose
adjustable **subcomponents** that enable more advanced workflows.

### Example: Sequential Environment Parameters

```julia
param = BendersSeqParam(
    iter_limit = 500,
    time_limit = 1800.0,
    verbose = true,
)

env = BendersSeq(master, oracle; param = param)
```

Typical environment parameters include:

* `iter_limit`: maximum number of iterations,
* `halt_limit`: maximum number of iterations without significant improvement,
* `time_limit`: global time limit,
* `gap_tolerance`: convergence tolerance,
* `verbose`: logging verbosity.

See [`BendersSeqInOutParam`](@ref) and [`BendersBnBParam`](@ref) for parameters
specific to `BendersSeqInOut` and `BendersBnB`, respectively.

---

## Benders Branch-and-Bound with Root Preprocessing and Callbacks

In addition to parameters defined in [`BendersBnBParam`](@ref),
[`BendersBnB`](@ref) supports configurable **root-node preprocessing** and
**callback-driven** workflows.

### 1) Simple lazy-callback–based usage (typical oracle, no user callback)

#### Using the default constructor

The short form `BendersBnB(master, oracle; param = benders_param)` wires the
oracle as a lazy callback and uses no root preprocessing:

```julia
master = Master(data; customize = customize_master_model!)
oracle = ClassicalOracle(data, master; customize = customize_sub_model!)

env = BendersBnB(master, oracle; param = benders_param)
log = solve!(env)
```

#### Explicit form (NoSeq / Seq / SeqInOut)

When finer control is needed, the environment can be constructed explicitly:

```julia
# NoSeq: no root Benders run
root_preprocessing = NoRootNodePreprocessing()
lazy_callback = LazyCallback(oracle)
user_callback = NoUserCallback()

env = BendersBnB(master, root_preprocessing, lazy_callback, user_callback; param = benders_param)

# Seq: sequential Benders at the root node
root_preprocessing = RootNodePreprocessing(
    oracle, BendersSeq,
    BendersSeqParam(time_limit = 200.0, gap_tolerance = 1e-9),
)

env = BendersBnB(master, root_preprocessing, LazyCallback(oracle), NoUserCallback(); param = benders_param)

# SeqInOut: stabilized In–Out Benders at the root node
root_preprocessing = RootNodePreprocessing(
    oracle, BendersSeqInOut,
    BendersSeqInOutParam(
        time_limit = 300.0,
        gap_tolerance = 1e-9,
        stabilizing_x = ones(data.n_facilities),
        α = 0.9,
        λ = 0.1,
    ),
)

env = BendersBnB(master, root_preprocessing, LazyCallback(oracle), NoUserCallback(); param = benders_param)
```

**When to use**

* `NoSeq`: no need for root tightening.
* `Seq`: improved root bound via sequential Benders.
* `SeqInOut`: stabilized root solve when a good `stabilizing_x` is available.

---

### 2) User-callback–based usage (e.g., disjunctive separation)

When stronger cuts are needed—such as separating fractional solutions—users can
add a **user callback**, for example with a disjunctive oracle.

```julia
using CPLEX

# typical oracles (κ, ν)
kappa = ClassicalOracle(data, master; customize = customize_sub_model!)
nu    = ClassicalOracle(data, master; customize = customize_sub_model!)
typical_oracles = [kappa, nu]

# DCGLP and SplitOracle parameters
dcglp_optimizer = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => 7,
    MOI.Silent() => true,
)
dcglp_param = DcglpParam(dcglp_optimizer; time_limit = 200.0)
split_param = SplitOracleParam(
    dcglp_param;
    norm = LpNorm(1.0),
    strengthened = true,
    lift = true,
    add_benders_cuts_to_master = 1,
)

# disjunctive oracle and callbacks
disjunctive_oracle = SplitOracle(master, typical_oracles, split_param)
user_callback = UserCallback(disjunctive_oracle; params = UserCallbackParam(frequency = 1))

lazy_oracle = ClassicalOracle(data, master; customize = customize_sub_model!)
lazy_callback = LazyCallback(lazy_oracle)

env = BendersBnB(
    master,
    NoRootNodePreprocessing(),
    lazy_callback,
    user_callback;
    param = benders_param,
)
log = solve!(env)
```

The master and DCGLP optimizers should be attached through standard JuMP APIs.

---

## Choosing the Right Environment

| Environment type            | Typical use case                     |
| --------------------------- | ------------------------------------ |
| `BendersSeq`                | Standard sequential Benders          |
| `BendersSeqInOut`           | In–Out stabilized sequential Benders |
| `BendersBnB`                | MILP master with callbacks           |
| Custom `AbstractBendersEnv` | Research and algorithm prototyping   |

---

## Summary

* Environments define *how* Benders decomposition is executed
* Swapping environments requires minimal code changes
* Behavior is controlled via environment-specific parameters and subcomponents
* Advanced workflows are enabled by combining root preprocessing and callbacks

This separation between **modeling**, **oracles**, and **execution environments**
allows BendersX to support a wide range of algorithms while maintaining a clean,
extensible API.
