include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV


#-----------------------------------------------------------------------
solver = :Gurobi


settings = SplitBenders.parse_commandline()
instance = "f200-c200-r5.0-p2"
data = SplitBenders.read_random_data(instance)

# instance = "p35"
# data = SplitBenders.read_data(instance)

#-----------------------------------------------------------------------
algo_params = SplitBenders.AlgorithmParams()

# "SPLIT_CUTSTRATEGY"
cut_strategy = "SPLIT_CUTSTRATEGY"

# "L1GAMMANORM", "L2GAMMANORM", "LINFGAMMANORM" "STANDARDNORM"
SplitCGLPNormType = "L1GAMMANORM"

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


master_env = SplitBenders.MasterProblem(data)
relax_integrality(master_env.model)
sub_env = SplitBenders.CFLPSplitSubEnv(data,algo_params)

df = SplitBenders.run_Benders(data,master_env,sub_env)

# result post processing
CSV.write("temp/result_$(instance)_$(cut_strategy)_$(SplitCGLPNormType)_$(SplitSetSelectionPolicy)_$(StrengthenCutStrategy)_$(SplitBendersStrategy).csv", df)





#--------------------------------------------------------------
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
# sub_env = SplitBenders.CFLPSplitSubEnv(data,algo_params)

# df = SplitBenders.run_Benders(data,master_env,sub_env)

# # result post processing
# cut_strategy = SplitBenders.reConvert_cutstrategy(cut_strategy)
# SplitCGLPNormType = SplitBenders.reConvert_normtype(SplitCGLPNormType)
# SplitSetSelectionPolicy = SplitBenders.reConvert_splitsetselectionpolicy(SplitSetSelectionPolicy)
# StrengthenCutStrategy = SplitBenders.reConvert_splitstengtheningpolicy(StrengthenCutStrategy)
# SplitBendersStrategy = SplitBenders.reConvert_splitbenderspolicy(SplitBendersStrategy)
# CSV.write("results/Gurobi/result_$(instance)_$(cut_strategy)_$(SplitCGLPNormType)_$(SplitSetSelectionPolicy)_$(StrengthenCutStrategy)_$(SplitBendersStrategy).csv", df)
