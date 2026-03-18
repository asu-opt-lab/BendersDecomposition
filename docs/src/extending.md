```@meta
CurrentModule = BendersX
```

# [Extending BendersX](@id extending-bendersx)

This page documents the public extension points that are meant to be reused in
custom algorithm development.

## Custom modeling hooks

The lowest-friction extension path is to reuse the built-in algorithms with
your own modeling code:

- define a subtype of [`AbstractData`](@ref);
- implement [`customize_master_model!`](@ref);
- implement [`customize_sub_model!`](@ref);
- optionally implement [`customize_mip_model!`](@ref) for direct-baseline
  comparisons.

This is the right starting point whenever the algorithmic workflow stays the
same and only the underlying optimization model changes.

## Custom oracles

Create a custom oracle when you want to change **which cuts are generated** or
how subproblems are evaluated.

```julia
import BendersX: AbstractOracle, generate_cuts

struct MyOracle <: AbstractOracle
    # state, parameters, models, caches, ...
end

function generate_cuts(
    oracle::MyOracle,
    x_value::Vector{Float64},
    t_value::Vector{Float64};
    tol_normalize = 1.0,
    time_limit = 3600.0,
)
    # return (is_in_L, hyperplanes, f_x)
end
```

`generate_cuts` must return:

- `is_in_L::Bool`: whether the candidate point is already accepted;
- `hyperplanes::Vector{Hyperplane}`: cuts to add to the master;
- `f_x::Vector{Float64}`: recourse evaluations used for upper-bound updates.

Useful helpers:

- [`Hyperplane`](@ref)
- [`aggregate`](@ref)
- [`evaluate_violation`](@ref)
- [`select_top_fraction`](@ref)
- [`hyperplanes_to_expression`](@ref)
- [`set_parameter!`](@ref)

## Custom environments

Create a custom environment when you want to change **how the algorithm runs**.

```julia
import BendersX: AbstractBendersEnv, solve!

struct MyEnv <: AbstractBendersEnv
    master::Master
    oracle::AbstractOracle
end

function solve!(env::MyEnv)
    # custom control flow
end
```

Custom environments are the right place for:

- alternative solve loops;
- hybrid sequential/callback workflows;
- custom stabilization;
- experiment-specific logging or termination logic.

## Generalized bound constraints

Typical oracles support generalized bound constraints through the return value
of [`customize_sub_model!`](@ref). The contract is:

```julia
return gbc_lhs, gbc_rhs, gbc_sense
```

where:

- `gbc_lhs` contains subproblem variables;
- `gbc_rhs` contains copied master variables or affine expressions in those
  variables;
- `gbc_sense` contains `UpperBound`, `LowerBound`, or `Fixed`.

This is the mechanism to express coupling relations that should be re-evaluated
whenever the oracle changes the fixed master point.

## Structural helpers

These helpers are public because they are useful when writing custom oracles or
specialized model transformations:

- [`copy_variables!`](@ref) mirrors master variable containers inside another
  JuMP model.
- [`var_from_tuple`](@ref) flattens a `NamedTuple` of JuMP containers into a
  `Vector{VariableRef}`.
- [`add_constraints`](@ref) appends batches of cut expressions under a symbolic
  model key.
- [`transfer_scaled_linear_rows_and_bounds_with_types!`](@ref) copies eligible
  master rows into the DCGLP formulation used by split cuts.

## Debugging and diagnostics

- [`infeasibility_report`](@ref) helps verify candidate master solutions.
- [`assign_attributes!`](@ref) applies a dictionary of solver attributes to a
  JuMP model.

When extending the package, start with a sequential environment and only add
callback or disjunctive machinery once the model-level contracts are stable.
