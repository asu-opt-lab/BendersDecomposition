```@meta
CurrentModule = BendersX
```

# [API Reference](@id api)

This page groups the documented public surface of `BendersX`. It includes both
exported names and public-but-not-auto-imported extension points.

## Core workflow

```@docs
AbstractData
AbstractMaster
Master
customize_master_model!
customize_sub_model!
customize_mip_model!
solve!
```

## Oracles

```@docs
AbstractOracle
AbstractTypicalOracle
AbstractOracleParam
BasicOracleParam
ClassicalOracleParam
ClassicalOracle
UnifiedOracleParam
UnifiedOracle
ParetoOracleParam
ParetoOracle
SeparableOracleParam
SeparableOracle
SplitOracleParam
DcglpParam
SplitOracle
CFLKnapsackOracleParam
CFLKnapsackOracle
UFLKnapsackOracleParam
UFLKnapsackOracle
generate_cuts
set_parameter!
```

## Environments

```@docs
AbstractBendersEnv
AbstractBendersSeq
AbstractBendersBnB
BendersSeqParam
BendersSeq
BendersSeqInOutParam
BendersSeqInOut
BendersBnBParam
BendersBnB
SpecializedBendersSeqParam
SpecializedBendersSeq
```

## Callbacks and preprocessing

```@docs
AbstractRootNodePreprocessing
NoRootNodePreprocessing
RootNodePreprocessing
DisjunctiveRootNodePreprocessing
LazyCallback
UserCallbackParam
NoUserCallback
UserCallback
```

## Problem helpers

```@docs
CFLPData
UFLPData
SCFLPData
SNIPData
read_GK_data
read_cfl_file
read_cflp_benchmark_data
read_uflp_benchmark_data
read_Simple_data
read_stochastic_capacited_facility_location_problem
read_snip_data
```

## Utilities

```@docs
GBCBoundType
TerminationStatus
NotSolved
TimeLimit
Optimal
InfeasibleOrNumericalIssue
LpNorm
SplitIndexSelectionRule
RandomFractional
MostFractional
LargestFractional
DisjunctiveCutsAppendRule
NoDisjunctiveCuts
AllDisjunctiveCuts
DisjunctiveCutsSmallerIndices
Hyperplane
aggregate
evaluate_violation
select_top_fraction
hyperplanes_to_expression
add_constraints
copy_variables!
var_from_tuple
assign_attributes!
transfer_scaled_linear_rows_and_bounds_with_types!
infeasibility_report
```
