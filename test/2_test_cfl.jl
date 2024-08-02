include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV, Logging


#-----------------------------------------------------------------------
# solver = :Gurobi
solver = :CPLEX

settings = SplitBenders.parse_commandline()
# instance = settings["instance"]
# instance = "f700-c700-r5.0-p3"
# data = SplitBenders.read_GK_data(instance)

instance = "p71"
data = SplitBenders.read_benchmark_data(instance)

#-----------------------------------------------------------------------
algo_params = SplitBenders.AlgorithmParams()

# "SPLIT_CUTSTRATEGY"
cut_strategy = "SPLIT_CUTSTRATEGY"

# "L1GAMMANORM", "L2GAMMANORM", "LINFGAMMANORM" "STANDARDNORM"
SplitCGLPNormType = "STANDARDNORM"

# "MOST_FRAC_INDEX", "RANDOM_INDEX"
SplitSetSelectionPolicy = "MOST_FRAC_INDEX"

# "SPLIT_PURE_CUT_STRATEGY", "SPLIT_STRENGTHEN_CUT_STRATEGY"
StrengthenCutStrategy = "SPLIT_STRENGTHEN_CUT_STRATEGY"

# "NO_SPLIT_BENDERS_STRATEGY", "ALL_SPLIT_BENDERS_STRATEGY", "TIGHT_SPLIT_BENDERS_STRATEGY"
SplitBendersStrategy = "ALL_SPLIT_BENDERS_STRATEGY"


SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractCutStrategy, cut_strategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractNormType, SplitCGLPNormType)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitSetSelectionPolicy, SplitSetSelectionPolicy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitStengtheningPolicy, StrengthenCutStrategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitBendersPolicy, SplitBendersStrategy)

master_env = SplitBenders.CFLPMasterProblem(data, solver=solver)
relax_integrality(master_env.model)
sub_env = SplitBenders.CFLPSplitSubEnv(data,algo_params, solver=solver)

df = SplitBenders.run_Benders(data,master_env,sub_env)





