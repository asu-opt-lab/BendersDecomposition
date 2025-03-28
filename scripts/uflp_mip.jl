using JuMP, DataFrames,Logging, CSV, CPLEX
using BendersDecomposition
import Random

Random.seed!(1218)

instance, output_dir, cut_strategy, benders_params = load_all_you_need()
data = read_Simple_data(instance)
mip_env = create_milp(data)
# relax_integrality(mip_env.model)

set_optimizer(mip_env.model, CPLEX.Optimizer)
set_time_limit_sec(mip_env.model, 13600.0)
set_optimizer_attribute(mip_env.model, "CPX_PARAM_EPINT", 0.0)
set_optimizer_attribute(mip_env.model, "CPX_PARAM_EPGAP", 1e-9)
set_optimizer_attribute(mip_env.model, "CPX_PARAM_EPRHS", 1e-9)
MOI.set(mip_env.model, MOI.RelativeGapTolerance(), 1e-9) 
set_optimizer_attribute(mip_env.model, "CPXPARAM_MIP_Display", 3)

start_time = time()
set_optimizer_attribute(mip_env.model, MOI.Silent(), false)
JuMP.optimize!(mip_env.model)
solution_summary(mip_env.model)

# df_mip = DataFrame(
#     objective_value = JuMP.objective_value(mip_env.model),
#     termination_status = termination_status(mip_env.model)
# )

df_mip = DataFrame(
    elapsed_time = time() - start_time,
    node_count = JuMP.node_count(mip_env.model),
    objective_bound = JuMP.objective_bound(mip_env.model),
    objective_value = JuMP.objective_value(mip_env.model),
    termination_status = termination_status(mip_env.model)
)

CSV.write(output_dir * "/$(instance)_df_mip.csv", df_mip)