using JuMP
using CPLEX

sub_optimizer = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => 7,
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPOPT" => 1e-9,
    "CPX_PARAM_NUMERICALEMPHASIS" => 1,
    MOI.Silent() => true,
)

master_optimizer = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => 7,
    "CPX_PARAM_EPINT" => 1e-9,
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPGAP" => 1e-6,
    MOI.Silent() => true,
)

mip_optimizer = optimizer_with_attributes(
    CPLEX.Optimizer,
    "CPXPARAM_Threads" => 7,
    "CPX_PARAM_EPINT" => 1e-9,
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPGAP" => 1e-6,
    MOI.Silent() => true,
)
