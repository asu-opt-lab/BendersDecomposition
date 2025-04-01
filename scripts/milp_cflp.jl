using JuMP, DataFrames, Logging, CSV, ArgParse, CPLEX
using BendersDecomposition


# Parse command line arguments
instance, output_dir, cut_strategy, benders_params = load_all_you_need()
data = read_GK_data(instance)
milp = create_milp(data)
# relax_integrality.(milp.model)
set_optimizer_attribute(milp.model, "CPX_PARAM_EPINT", 0.0)
set_optimizer_attribute(milp.model, "CPX_PARAM_EPGAP", 1e-9)
set_optimizer_attribute(milp.model, "CPX_PARAM_EPRHS", 1e-9)
set_optimizer(milp.model, CPLEX.Optimizer)
set_time_limit_sec(milp.model, 3600.0)
optimize!(milp.model)
@info termination_status(milp.model)
@info "Solve time: $(solve_time(milp.model))"
@info "Objective value: $(objective_value(milp.model))"
@info "objective bound: $(objective_bound(milp.model))"
# @info "relative gap: $(relative_gap(milp.model))"
