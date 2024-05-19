include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV

settings = SplitBenders.parse_commandline()
instance = settings["instance"]
data = SplitBenders.read_random_data(instance)

# instance = "p2"
# data = SplitBenders.read_data(instance)

cut_strategy = settings["cut_strategy"]
SplitCGLPNormType = settings["SplitCGLPNormType"]
SplitSetSelectionPolicy = settings["SplitSetSelectionPolicy"]
StrengthenCutStrategy = settings["StrengthenCutStrategy"]
SplitBendersStrategy = settings["SplitBendersStrategy"]

algo_params = SplitBenders.AlgorithmParams(
    cut_strategy = c_cut_strategy,
    SplitCGLPNormType = c_SplitCGLPNormType,
    SplitSetSelectionPolicy = c_SplitSetSelectionPolicy,
    StrengthenCutStrategy = c_StrengthenCutStrategy,
    SplitBendersStrategy = c_SplitBendersStrategy
)

master_env = SplitBenders.MasterProblem(data)
sub_env = SplitBenders.CFLPStandardSubEnv(data,algo_params)

df = SplitBenders.run_Benders(data,master_env,sub_env)

# result post processing
cut_strategy = SplitBenders.reConvert_cutstrategy(cut_strategy)
SplitCGLPNormType = SplitBenders.reConvert_normtype(SplitCGLPNormType)
SplitSetSelectionPolicy = SplitBenders.reConvert_splitsetselectionpolicy(SplitSetSelectionPolicy)
StrengthenCutStrategy = SplitBenders.reConvert_splitstengtheningpolicy(StrengthenCutStrategy)
SplitBendersStrategy = SplitBenders.reConvert_splitbenderspolicy(SplitBendersStrategy)
CSV.write("results/Gurobi/result_$(instance)_$(cut_strategy)_$(SplitCGLPNormType)_$(SplitSetSelectionPolicy)_$(StrengthenCutStrategy)_$(SplitBendersStrategy).csv", df)
