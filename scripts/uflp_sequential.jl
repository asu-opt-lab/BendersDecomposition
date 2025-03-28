
using JuMP, DataFrames,Logging, CSV
using BendersDecomposition
import Random

Random.seed!(1218)

instance, output_dir, cut_strategy, benders_params = load_all_you_need()
data = read_Simple_data(instance)
loop_strategy = Sequential()
@info "Running Benders decomposition"
df_root_node_preprocessing = run_Benders(data, loop_strategy, cut_strategy, benders_params)

@info "Writing results to CSV"
CSV.write(output_dir * "/$(instance)_df_root_node_preprocessing.csv", df_root_node_preprocessing)
