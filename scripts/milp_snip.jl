using JuMP, DataFrames, Logging, CSV, ArgParse, CPLEX
using BendersDecomposition


# Parse command line arguments
instance, snip_no, budget, output_dir, cut_strategy, benders_params = load_snip_data()

# Load other necessary components
data = read_snip_data(instance, snip_no, budget)
milp = create_milp(data)
# relax_integrality.(milp.model)
set_optimizer(milp.model, CPLEX.Optimizer)
set_time_limit_sec(milp.model, 3600.0)
optimize!(milp.model)
@info termination_status(milp.model)
@info "Solve time: $(solve_time(milp.model))"
@info "Objective value: $(objective_value(milp.model))"

