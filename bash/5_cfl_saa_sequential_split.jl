include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV

settings = SplitBenders.parse_commandline()
nscenario = 10
datas = Array{SplitBenders.CFLPData}(undef, nscenario)
for w in 1:nscenario
    instance = "saa-f500-c500-r5.0-p$w"
    datas[w] = SplitBenders.read_random_data(instance)
end


# instance = "p70"
# data = SplitBenders.read_data(instance)

#-----------------------------------------------------------------------
algo_params = SplitBenders.AlgorithmParams()

# "SPLIT_CUTSTRATEGY"
cut_strategy = "SPLIT_CUTSTRATEGY"

# "L1GAMMANORM", "L2GAMMANORM", "LINFGAMMANORM" "STANDARDNORM"
SplitCGLPNormType = "LINFGAMMANORM"

# "MOST_FRAC_INDEX", "RANDOM_INDEX"
SplitSetSelectionPolicy = "RANDOM_INDEX"

# "SPLIT_PURE_CUT_STRATEGY", "SPLIT_STRENGTHEN_CUT_STRATEGY"
StrengthenCutStrategy = "SPLIT_STRENGTHEN_CUT_STRATEGY"

# "NO_SPLIT_BENDERS_STRATEGY", "ALL_SPLIT_BENDERS_STRATEGY", "TIGHT_SPLIT_BENDERS_STRATEGY"
SplitBendersStrategy = "TIGHT_SPLIT_BENDERS_STRATEGY"


SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractCutStrategy, cut_strategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractNormType, SplitCGLPNormType)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitSetSelectionPolicy, SplitSetSelectionPolicy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitStengtheningPolicy, StrengthenCutStrategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitBendersPolicy, SplitBendersStrategy)


master_env = SplitBenders.MasterProblem(datas,nscenario)
undo = relax_integrality(master_env.model)
sub_envs = Array{SplitBenders.CFLPSplitSubEnv}(undef, nscenario)
for w in 1:nscenario
    sub_envs[w] = SplitBenders.CFLPSplitSubEnv(datas[w],algo_params)
end

df = SplitBenders.run_Benders_SAA(datas,master_env,sub_envs,nscenario)

undo()

optimize!(master_env.model)

@info objective_value(master_env.model)