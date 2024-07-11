include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV

# solver = :Gurobi
solver = :CPLEX
time_limit = 1000

settings = SplitBenders.parse_commandline()
# instance = settings["instance"]
# data = SplitBenders.read_random_data(instance)
# instance = "f700-c700-r5.0-p10"
# data = SplitBenders.read_random_data(instance)

instance = "p70"
data = SplitBenders.read_data(instance)

#-----------------------------------------------------------------------
algo_params = SplitBenders.AlgorithmParams()

# "ORDINARY_CUTSTRATEGY" "ADVANCED_CUTSTRATEGY" "KN_CUTSTRATEGY"
cut_strategy = "ADVANCED_CUTSTRATEGY"

# "L1GAMMANORM", "L2GAMMANORM", "LINFGAMMANORM" "STANDARDNORM"
SplitCGLPNormType = "nothing"

# "MOST_FRAC_INDEX", "RANDOM_INDEX"
SplitSetSelectionPolicy = "nothing"

# "SPLIT_PURE_CUT_STRATEGY", "SPLIT_STRENGTHEN_CUT_STRATEGY"
StrengthenCutStrategy = "nothing"

# "NO_SPLIT_BENDERS_STRATEGY", "ALL_SPLIT_BENDERS_STRATEGY", "TIGHT_SPLIT_BENDERS_STRATEGY"
SplitBendersStrategy = "nothing"


SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractCutStrategy, cut_strategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractNormType, SplitCGLPNormType)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitSetSelectionPolicy, SplitSetSelectionPolicy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitStengtheningPolicy, StrengthenCutStrategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitBendersPolicy, SplitBendersStrategy)


master_env = SplitBenders.MasterProblem(data, solver = solver)
relax_integrality(master_env.model)
# sub_env = SplitBenders.CFLPStandardSubEnv(data,algo_params, solver = solver)
sub_env = SplitBenders.CFLPStandardADSubEnv(data,algo_params, solver = solver)
# sub_env = SplitBenders.CFLPStandardKNSubEnv(data,algo_params, solver = solver)

df = SplitBenders.run_Benders(data,master_env,sub_env)

# result post processing
# CSV.write("results/Gurobi/result_$(instance)_$(cut_strategy)_$(SplitCGLPNormType)_$(SplitSetSelectionPolicy)_$(StrengthenCutStrategy)_$(SplitBendersStrategy).csv", df)
# CSV.write("results2/Ordinary/result_$(instance).csv", df)
# CSV.write("results2/Advanced/result_$(instance).csv", df)


# cut_strategy = settings["cut_strategy"]
# SplitCGLPNormType = settings["SplitCGLPNormType"]
# SplitSetSelectionPolicy = settings["SplitSetSelectionPolicy"]
# StrengthenCutStrategy = settings["StrengthenCutStrategy"]
# SplitBendersStrategy = settings["SplitBendersStrategy"]

# algo_params = SplitBenders.AlgorithmParams(
#     cut_strategy = cut_strategy,
#     SplitCGLPNormType = SplitCGLPNormType,
#     SplitSetSelectionPolicy = SplitSetSelectionPolicy,
#     StrengthenCutStrategy = StrengthenCutStrategy,
#     SplitBendersStrategy = SplitBendersStrategy
# )

# master_env = SplitBenders.MasterProblem(data)
# sub_env = SplitBenders.CFLPStandardSubEnv(data,algo_params)

# df = SplitBenders.run_Benders(data,master_env,sub_env)

# # result post processing
# cut_strategy = SplitBenders.reConvert_cutstrategy(cut_strategy)
# SplitCGLPNormType = SplitBenders.reConvert_normtype(SplitCGLPNormType)
# SplitSetSelectionPolicy = SplitBenders.reConvert_splitsetselectionpolicy(SplitSetSelectionPolicy)
# StrengthenCutStrategy = SplitBenders.reConvert_splitstengtheningpolicy(StrengthenCutStrategy)
# SplitBendersStrategy = SplitBenders.reConvert_splitbenderspolicy(SplitBendersStrategy)
# CSV.write("results/Gurobi/result_$(instance)_$(cut_strategy)_$(SplitCGLPNormType)_$(SplitSetSelectionPolicy)_$(StrengthenCutStrategy)_$(SplitBendersStrategy).csv", df)
