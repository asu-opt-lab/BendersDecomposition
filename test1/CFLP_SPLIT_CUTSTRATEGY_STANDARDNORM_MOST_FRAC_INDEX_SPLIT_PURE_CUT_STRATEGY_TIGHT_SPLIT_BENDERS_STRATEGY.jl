include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV

settings = SplitBenders.parse_commandline()
instance = settings["instance"]
data = SplitBenders.read_random_data(instance)

algo_params = SplitBenders.AlgorithmParams(
    cut_strategy=Main.SplitBenders.SplitCutStrategy(), 
    SplitCGLPNormType=Main.SplitBenders.StandardNorm(), 
    SplitSetSelectionPolicy=Main.SplitBenders.MostFracIndex(), 
    StrengthenCutStrategy=Main.SplitBenders.SplitPureCutStrategy(), 
    SplitBendersStrategy=Main.SplitBenders.TightSplitBendersStrategy()
)

master_env = SplitBenders.MasterProblem(data)
sub_env = SplitBenders.CFLPSplitSubEnv(data,algo_params)

df = SplitBenders.run_Benders(data,master_env,sub_env)

CSV.write("results/Gurobi/result_$(instance)_SPLIT_CUTSTRATEGY_STANDARDNORM_MOST_FRAC_INDEX_SPLIT_PURE_CUT_STRATEGY_TIGHT_SPLIT_BENDERS_STRATEGY.csv", df)
