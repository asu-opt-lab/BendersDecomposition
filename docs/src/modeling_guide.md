```@meta
CurrentModule = BendersX
```

# [Modeling Guide](@id modeling-guide)

BendersX keeps the algorithmic engine separate from the model definition. You
provide the master and subproblem formulations through standard JuMP code.

## Data containers

Every modeling workflow starts with a subtype of [`AbstractData`](@ref). The
package does not constrain the field layout; the data object only needs to
carry the information required to build your master and subproblem models.

```julia
struct MyData <: AbstractData
    # instance-specific fields
end
```

## Master contract

Implement [`customize_master_model!`](@ref) for your data type:

```julia
customize_master_model!(model::Model, data::MyData) -> NamedTuple, Vector{VariableRef}
```

The function must:

1. choose and configure the optimizer;
2. add the first-stage variables and constraints;
3. define the master objective;
4. return the coupling variables as a `NamedTuple`;
5. return the approximation variable or variables `t`.

```julia
function customize_master_model!(model::Model, data::MyData)
    # build x and t
    return (x = x, u = u), t
end
```

The names in the returned `NamedTuple` matter: they become the keyword names
received by [`customize_sub_model!`](@ref).

## Subproblem contract

Implement [`customize_sub_model!`](@ref) for the same data type:

```julia
customize_sub_model!(model::Model, data::MyData, scen_idx::Int; kwargs...)
```

The function must:

- build the subproblem variables and constraints;
- use the copied master variables passed through keyword arguments;
- return `nothing` for the standard case, or a generalized bound constraint
  triple for advanced use.

```julia
function customize_sub_model!(model::Model, data::MyData, scen_idx::Int; x, u)
    @variable(model, y >= 0)
    @constraint(model, y <= x[1] + u[2])
    return nothing
end
```

!!! warning "Do not redeclare master variables"
    The keyword arguments passed into `customize_sub_model!` are already the
    copied master variables inside the subproblem model. Reuse them directly;
    do not add a second set of coupling variables.

## Variable containers

`BendersX` supports scalar variables together with the common JuMP containers:

- plain arrays;
- `DenseAxisArray`;
- `SparseAxisArray`.

[`copy_variables!`](@ref) and [`var_from_tuple`](@ref) preserve that structure
when the oracle copies master variables into the subproblem.

## Generalized bound constraints

Typical oracles support an optional generalized bound constraint contract. In
that case, `customize_sub_model!` may return:

```julia
(gbc_lhs, gbc_rhs, gbc_sense)
```

with:

- `gbc_lhs::Vector{VariableRef}` containing **subproblem** variables only;
- `gbc_rhs::Vector{Union{VariableRef,AffExpr}}` containing affine expressions in
  copied master variables only;
- `gbc_sense::Vector{GBCBoundType}` containing `UpperBound`, `LowerBound`, or
  `Fixed`.

Example:

```julia
function customize_sub_model!(model::Model, data::MyData, scen_idx::Int; x)
    @variable(model, y >= 0)
    @constraint(model, y >= 0)

    gbc_lhs = [y]
    gbc_rhs = [2.0 * x[1] + x[2]]
    gbc_sense = [UpperBound]

    return gbc_lhs, gbc_rhs, gbc_sense
end
```

Use this only when the oracle should enforce bounds of the form
`subproblem_var <= affine(master_vars)`, `>=`, or `==` while evaluating a
candidate point.

## Direct MIP baselines

[`customize_mip_model!`](@ref) is the public hook for building a direct
monolithic MIP baseline. The built-in problem libraries provide methods for the
shipped data types, which is useful when you want to compare a Benders solve
against a direct solve of the full formulation.

## Practical checklist

- Make the keyword names in [`customize_sub_model!`](@ref) match the names
  returned by [`customize_master_model!`](@ref).
- Keep typical subproblems LP-compatible if you plan to use
  [`ClassicalOracle`](@ref), [`UnifiedOracle`](@ref), or [`ParetoOracle`](@ref).
- Use `scen_idx` even when the first version of your model is not
  scenario-based; it keeps the interface compatible with separable workflows.
- Start with a sequential environment and only move to callback-based
  environments after the model-level contract is behaving as expected.
