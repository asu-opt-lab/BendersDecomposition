# PolarDCGLP Infeasibility Analysis

## Context

The infeasibility is in the auxiliary `PolarDCGLP` model, not in the original master problem itself. The IIS indicates that the conflict is mainly caused by the interaction among:

- `conx`
- `coneta`
- `condelta`
- the split constraints `con_split_kappa` and `con_split_nu`
- one transferred master constraint, which appears to be the global capacity constraint

This means the issue is structural: the current split is incompatible with the scaled block-wise feasibility requirements of the polar model at the current point.

## Key mechanism

For the current IIS diagnostic:

- `j_split = 183`
- `phi_0 = 0`
- `x_value[183] = 0.09375`

The main polar constraints are created in `polardcglp/src/PolarDCGLP.jl`:

```julia
@constraint(dcglp, coneta[i in 1:2, j in 1:master.dim_x], 0 >= -omega_0[i] + omega_x[i, j])
@constraint(dcglp, condelta[i in 1:2, j in 1:master.dim_x], 0 >= -omega_x[i, j])
@constraint(dcglp, con0, omega_0[1] + omega_0[2] == 1.0)
@constraint(dcglp, conx[j in 1:master.dim_x], omega_x[1, j] + omega_x[2, j] == 0.0)
```

Then in `polardcglp/src/PolarDCGLPInterface.jl`, the right-hand side of `conx` is replaced by the current point:

```julia
JuMP.set_normalized_rhs.(oracle.dcglp[:conx], x_value)
```

So the actual constraint used during the solve is:

```text
conx[j]: omega_x[1,j] + omega_x[2,j] = x_value[j]
```

Also, `coneta` and `condelta` mean:

```text
coneta[i,j]:   omega_x[i,j] <= omega_0[i]
condelta[i,j]: omega_x[i,j] >= 0
```

The split constraints are created in `polardcglp/src/PolarDCGLP.jl`:

```julia
@constraint(dcglp, con_split_kappa, 0 >= dcglp[:omega_0][1] * (phi_0 + 1.0) - phi' * dcglp[:omega_x][1, :])
@constraint(dcglp, con_split_nu, 0 >= -dcglp[:omega_0][2] * phi_0 + phi' * dcglp[:omega_x][2, :])
```

For the selected split `phi = e_183`, `phi_0 = 0`, these become:

```text
con_split_kappa: 0 >= omega_0[1] - omega_x[1,183]
con_split_nu:    0 >= omega_x[2,183]
```

Now the exact derivation is:

```text
1. From condelta[2,183]:
   omega_x[2,183] >= 0

2. From con_split_nu:
   omega_x[2,183] <= 0

3. Therefore:
   omega_x[2,183] = 0

4. From conx[183]:
   omega_x[1,183] + omega_x[2,183] = x_value[183] = 0.09375

5. Therefore:
   omega_x[1,183] = 0.09375

6. From con_split_kappa:
   omega_x[1,183] >= omega_0[1]
   so omega_0[1] <= 0.09375

7. From coneta[1,183]:
   omega_x[1,183] <= omega_0[1]
   so omega_0[1] >= 0.09375

8. Hence:
   omega_0[1] = 0.09375
   omega_0[2] = 1 - omega_0[1] = 0.90625
```

So the split does not merely imply an inequality. It fixes the block weights exactly:

```text
omega_0[1] = x_value[183]
omega_0[2] = 1 - x_value[183]
```

In this case, the second polar block is forced to carry exactly `90.625%` of the total scaling weight.

## Why this becomes infeasible

The polar construction transfers each master linear constraint to each block in scaled form:

```text
A * omega_x[i,:] >= b * omega_0[i]
```

The transfer code in `src/utils/utilsInterface.jl` is:

```julia
if S == MOI.GreaterThan{Float64}
    @constraint(dcglp, expr + constant * omega0 >= set.lower * omega0)
elseif S == MOI.LessThan{Float64}
    @constraint(dcglp, expr + constant * omega0 <= set.upper * omega0)
elseif S == MOI.EqualTo{Float64}
    @constraint(dcglp, expr + constant * omega0 == set.value * omega0)
elseif S == MOI.Interval{Float64}
    @constraint(dcglp, expr + constant * omega0 <= set.upper * omega0)
    @constraint(dcglp, expr + constant * omega0 >= set.lower * omega0)
end
```

The anonymous IIS inequality

```text
-6397 omega_0[2] + 73 omega_x[2,1] + ... + 59 omega_x[2,300] >= 0
```

looks exactly like the master total-capacity constraint after scaling into block 2. In other words, block 2 must satisfy something of the form:

```text
sum(capacity_j * omega_x[2,j]) >= 6397 * omega_0[2]
```

For a CFLP instance, the original master constraint is:

```julia
@constraint(model, capacity, sum(data.capacities[i] * x[i] for i in 1:I) >= sum(data.demands))
```

For an SCFLP instance, the original master constraint is:

```julia
@constraint(model, capacity, sum(data.capacities[i] * x[i] for i in 1:I) >= max_demand)
```

But at the same time:

- each `omega_x[2,j]` is nonnegative,
- each `omega_x[2,j]` is limited by `conx[j]` and cannot exceed what is available from the current `x_value[j]`,
- and the split removes variable 183 from block 2 by forcing `omega_x[2,183] = 0`.

Therefore, block 2 is asked to satisfy a very large fraction of the total capacity requirement while being restricted to only a limited share of the current `x_value`. At the current point, that combination is impossible, so the polar model becomes infeasible.

## Interpretation

This does **not** mean that the original master problem is infeasible.

It also does **not** imply that a global valid cut such as `x[183] = 1` can be added directly.

What it means is:

- the current split on `x[183]` is incompatible with the polar block decomposition at this `x_value`,
- the split fixes a very unbalanced scaling (`omega_0[2] = 0.90625`),
- and the scaled master feasibility requirements for the second block can no longer be satisfied.

## Practical conclusion

The immediate cause of infeasibility is:

1. the split on `x[183]`,
2. the resulting exact weight `omega_0[2] = 0.90625`,
3. the transferred block-2 capacity requirement,
4. and the inability of `omega_x[2,:]` to supply enough scaled capacity under the current `x_value`.

So the failure is best interpreted as:

> this split is not viable for the current point in the polar DCGLP formulation.

