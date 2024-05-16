include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP

# instance = "f300-c300-r5.0-p1"
# data = SplitBenders.read_random_data(instance)

instance = "p35"
data = SplitBenders.read_data(instance)

algo_params = SplitBenders.AlgorithmParams(
    cut_strategy= SplitBenders.SPLIT_CUTSTRATEGY,
    SplitCGLPNormType=SplitBenders.LINFGAMMANORM,
    SplitSetSelectionPolicy=SplitBenders.RANDOM_INDEX,
    StrengthenCutStrategy=SplitBenders.SPLIT_STRENGTHEN_CUT_STRATEGY,
    SplitBendersStrategy=SplitBenders.ALL_SPLIT_BENDERS_STRATEGY
)


master_env = SplitBenders.MasterProblem(data)
sub_env = SplitBenders.CFLPSplitSubEnv(data,algo_params)
SplitBenders.run_Benders(data,master_env,sub_env)


# result post processing