include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP

# instance = "f100-c100-r5.0-p2"
# data = SplitBenders.read_random_data(instance)

instance = "p2"
data = SplitBenders.read_data(instance)

algo_params = SplitBenders.AlgorithmParams(
    SplitBenders.KN_CUTSTRATEGY,
    nothing,
    nothing,
    nothing,
    nothing
)

master_env = SplitBenders.MasterProblem(data)
sub_env = SplitBenders.CFLPStandardKNSubEnv(data,algo_params)

SplitBenders.run_Benders(data,master_env,sub_env)


# result post processing