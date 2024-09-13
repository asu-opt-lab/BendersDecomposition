include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV, Logging, DataFrames, CPLEX, Gurobi

solver = :Gurobi
# solver = :CPLEX

# instance = "MO1"
# @info "Instance: $instance"
# data = SplitBenders.read_Orlib_data(instance; filepath = "src/BendersDatasets/M/O")

instance = "ga500a-1"
data = SplitBenders.read_Simple_data(instance; filepath = "src/BendersDatasets/KoerkelGhosh-asym/")

# instance = "p10"
# data = SplitBenders.read_benchmark_data(instance)

#-----------------------------------------------------------------------
algo_params = SplitBenders.AlgorithmParams()
# "SPLIT_CUTSTRATEGY"
cut_strategy = "SPLIT_CUTSTRATEGY"
# "L1GAMMANORM", "L2GAMMANORM", "LINFGAMMANORM" "STANDARDNORM"
SplitCGLPNormType = "L2GAMMANORM"
# "MOST_FRAC_INDEX", "RANDOM_INDEX"
SplitSetSelectionPolicy = "MOST_FRAC_INDEX"
# "SPLIT_PURE_CUT_STRATEGY", "SPLIT_STRENGTHEN_CUT_STRATEGY"
StrengthenCutStrategy = "SPLIT_PURE_CUT_STRATEGY"
# "NO_SPLIT_BENDERS_STRATEGY", "ALL_SPLIT_BENDERS_STRATEGY", "TIGHT_SPLIT_BENDERS_STRATEGY"
SplitBendersStrategy = "NO_SPLIT_BENDERS_STRATEGY"
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractCutStrategy, cut_strategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractNormType, SplitCGLPNormType)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitSetSelectionPolicy, SplitSetSelectionPolicy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitStengtheningPolicy, StrengthenCutStrategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitBendersPolicy, SplitBendersStrategy)
#-----------------------------------------------------------------------

master_env = SplitBenders.UFLPMasterProblem(data, solver=solver)
relax_integrality(master_env.model)
sub_env = SplitBenders.UFLPSplitSubEnv(data,algo_params, solver=solver)
SplitBenders.run_Benders(data,master_env,sub_env)
# set_binary.(master_env.model[:x])
# SplitBenders.run_Benders_callback(data,master_env,sub_env)