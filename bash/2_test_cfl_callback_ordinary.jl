include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP

settings = SplitBenders.parse_commandline()
instance = "f200-c200-r5.0-p1"
data = SplitBenders.read_random_data(instance)

# instance = "p70"
# data = SplitBenders.read_data(instance)


algo_params = SplitBenders.AlgorithmParams()
cut_strategy = settings["cut_strategy"]
SplitCGLPNormType = settings["SplitCGLPNormType"]
SplitSetSelectionPolicy = settings["SplitSetSelectionPolicy"]
StrengthenCutStrategy = settings["StrengthenCutStrategy"]
SplitBendersStrategy = settings["SplitBendersStrategy"]
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractCutStrategy, cut_strategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractNormType, SplitCGLPNormType)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitSetSelectionPolicy, SplitSetSelectionPolicy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitStengtheningPolicy, StrengthenCutStrategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitBendersPolicy, SplitBendersStrategy)



master_env = SplitBenders.MasterProblem(data)
sub_env = SplitBenders.CFLPStandardSubEnv(data,algo_params)

SplitBenders.run_Benders_callback(data,master_env,sub_env)


# result post processing