include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV, CPLEX, Gurobi

s = "ORDINARY_CUTSTRATEGY"
# ss = eval(Meta.parse(s))

ss = SplitBenders.eval(Meta.parse(s))
@info ss
@info typeof(ss)
@info ss == SplitBenders.ORDINARY_CUTSTRATEGY
@info typeof(ss) == SplitBenders.AbstractCutStrategy

global algo_params = SplitBenders.AlgorithmParams()
@info algo_params.cut_strategy

SplitBenders.set_params_attribute(algo_params.cut_strategy, s)
# SplitBenders.set_params_attribute(algo_params, s)


@info algo_params.cut_strategy
@info algo_params.cut_strategy == SplitBenders.ORDINARY_CUTSTRATEGY

# instance = "p2"
# data = SplitBenders.read_data(instance)
# master_env = SplitBenders.MasterProblem(data)
# sub_env = SplitBenders.CFLPSplitSubEnv(data,algo_params)

# df = SplitBenders.run_Benders(data,master_env,sub_env)
