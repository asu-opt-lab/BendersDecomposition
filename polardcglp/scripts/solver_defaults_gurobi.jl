using JuMP
using Gurobi

sub_optimizer = optimizer_with_attributes(
    Gurobi.Optimizer,
    "Threads" => 7,
    "FeasibilityTol" => 1e-9,
    "OptimalityTol" => 1e-9,
    "NumericFocus" => 1,
    MOI.Silent() => true,
)

master_optimizer = optimizer_with_attributes(
    Gurobi.Optimizer,
    "Threads" => 7,
    "IntFeasTol" => 1e-9,
    "FeasibilityTol" => 1e-9,
    "MIPGap" => 1e-6,
    MOI.Silent() => true,
)

mip_optimizer = optimizer_with_attributes(
    Gurobi.Optimizer,
    "Threads" => 7,
    "IntFeasTol" => 1e-9,
    "FeasibilityTol" => 1e-9,
    "MIPGap" => 1e-6,
    MOI.Silent() => true,
)
