include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV, CPLEX, Gurobi

# instance = "f50-c50-r5.0-p1"
# data = DisjunctiveBenders.read_random_data(instance)

instance = "p30"
data = SplitBenders.read_data(instance)
@info instance

mip_env = SplitBenders.CFLPMipEnv(data)
optimize!(mip_env.model)
correct_answer = objective_value(mip_env.model)


for cut_strategy in [SplitBenders.ORDINARY_CUTSTRATEGY, SplitBenders.SPLIT_CUTSTRATEGY]
    for SplitCGLPNormType in [SplitBenders.STANDARDNORM, SplitBenders.L1GAMMANORM, SplitBenders.LINFGAMMANORM]
        for SplitSetSelectionPolicy in [SplitBenders.MOST_FRAC_INDEX, SplitBenders.RANDOM_INDEX]
            for StrengthenCutStrategy in [SplitBenders.SPLIT_PURE_CUT_STRATEGY, SplitBenders.SPLIT_STRENGTHEN_CUT_STRATEGY]
                for SplitBendersStrategy in [SplitBenders.NO_SPLIT_BENDERS_STRATEGY, SplitBenders.ALL_SPLIT_BENDERS_STRATEGY, SplitBenders.TIGHT_SPLIT_BENDERS_STRATEGY]
                    algo_params = SplitBenders.AlgorithmParams(
                        cut_strategy,
                        SplitCGLPNormType,
                        SplitSetSelectionPolicy,
                        StrengthenCutStrategy,
                        SplitBendersStrategy
                    )
                    @info algo_params
                    master_env = SplitBenders.MasterProblem(data)
                    sub_env = SplitBenders.CFLPSplitSubEnv(data,algo_params)
                    SplitBenders.run_Benders(data,master_env,sub_env)
                    @assert objective_value(master_env.model) ≈ correct_answer
                end
            end
        end
    end
    
end