include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV

settings = SplitBenders.parse_commandline()
instance = settings["instance"]
data = SplitBenders.read_random_data(instance)

algo_params = SplitBenders.AlgorithmParams(
    cut_strategy=Main.SplitBenders.SplitCutStrategy(), 
    SplitCGLPNormType=Main.SplitBenders.LInfGammaNorm(), 
    SplitSetSelectionPolicy=Main.SplitBenders.RandomIndex(), 
    StrengthenCutStrategy=Main.SplitBenders.SplitPureCutStrategy(), 
    SplitBendersStrategy=Main.SplitBenders.NoSplitBendersStrategy()
)

master_env = SplitBenders.MasterProblem(data)
sub_env = SplitBenders.CFLPSplitSubEnv(data,algo_params)

df = SplitBenders.run_Benders(data,master_env,sub_env)

CSV.write("results/Gurobi/result_$(instance)_SPLIT_CUTSTRATEGY_LINFGAMMANORM_RANDOM_INDEX_SPLIT_PURE_CUT_STRATEGY_NO_SPLIT_BENDERS_STRATEGY.csv", df)
