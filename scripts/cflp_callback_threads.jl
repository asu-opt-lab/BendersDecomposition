
using JuMP, DataFrames,Logging, CSV
using BendersDecomposition, Random

Random.seed!(1234)

instance, output_dir, cut_strategy, benders_params = load_all_you_need()
data = read_GK_data(instance)
loop_strategy = Callback()
@info "Running Benders decomposition"
# df_root_node_preprocessing, df_callback = run_Benders(data, loop_strategy, cut_strategy, benders_params)
run_Benders(data, loop_strategy, cut_strategy, benders_params)

# @info "Writing results to CSV"
# CSV.write(output_dir * "/$(instance)_df_root_node_preprocessing.csv", df_root_node_preprocessing)
# CSV.write(output_dir * "/$(instance)_df_callback.csv", df_callback)
