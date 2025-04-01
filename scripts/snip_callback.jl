
using JuMP, DataFrames,Logging, CSV
using BendersDecomposition

instance, snip_no, budget, output_dir, cut_strategy, benders_params = load_snip_data()
data = read_snip_data(instance, snip_no, budget)
loop_strategy = StochasticSequential()
@info "Running Benders decomposition"
# df_root_node_preprocessing, df_callback = run_Benders(data, loop_strategy, cut_strategy, benders_params)
run_Benders(data, loop_strategy, cut_strategy, benders_params)

# @info "Writing results to CSV"
# CSV.write(output_dir * "/$(instance)_df_root_node_preprocessing.csv", df_root_node_preprocessing)
# CSV.write(output_dir * "/$(instance)_df_callback.csv", df_callback)
