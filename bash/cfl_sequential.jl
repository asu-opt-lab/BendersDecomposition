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
    cut_strategy=cut_strategy, 
    SplitCGLPNormType=SplitCGLPNormType, 
    SplitSetSelectionPolicy=SplitSetSelectionPolicy, 
    StrengthenCutStrategy=StrengthenCutStrategy, 
    SplitBendersStrategy=SplitBendersStrategy
)

master_env = SplitBenders.MasterProblem(data)
sub_env = SplitBenders.CFLPStandardSubEnv(data,algo_params)

df = SplitBenders.run_Benders(data,master_env,sub_env)

# result post processing
CSV.write("result/Gurobi/result_$(instance)_$(cut_strategy)_$(SplitCGLPNormType)_$(SplitSetSelectionPolicy)_$(StrengthenCutStrategy)_$(SplitBendersStrategy).csv", df)
