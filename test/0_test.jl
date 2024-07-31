include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV, Logging, DataFrames, CPLEX, Gurobi


#-----------------------------------------------------------------------
# solver = :Gurobi
solver = :CPLEX

settings = SplitBenders.parse_commandline()

#-----------------------------------------------------------------------
algo_params = SplitBenders.AlgorithmParams()

# "SPLIT_CUTSTRATEGY"
cut_strategy = "SPLIT_CUTSTRATEGY"

# "L1GAMMANORM", "L2GAMMANORM", "LINFGAMMANORM" "STANDARDNORM"
SplitCGLPNormType = "L1GAMMANORM"

# "MOST_FRAC_INDEX", "RANDOM_INDEX"
SplitSetSelectionPolicy = "MOST_FRAC_INDEX"

# "SPLIT_PURE_CUT_STRATEGY", "SPLIT_STRENGTHEN_CUT_STRATEGY"
StrengthenCutStrategy = "SPLIT_PURE_CUT_STRATEGY"

# "NO_SPLIT_BENDERS_STRATEGY", "ALL_SPLIT_BENDERS_STRATEGY", "TIGHT_SPLIT_BENDERS_STRATEGY"
SplitBendersStrategy = "NO_SPLIT_BENDERS_STRATEGY"




SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractCutStrategy, cut_strategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractNormType, SplitCGLPNormType)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitSetSelectionPolicy, SplitSetSelectionPolicy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitStengtheningPolicy, StrengthenCutStrategy)
SplitBenders.set_params_attribute(algo_params, SplitBenders.AbstractSplitBendersPolicy, SplitBendersStrategy)

# master_env = SplitBenders.MasterProblem(data)
# relax_integrality(master_env.model)
# sub_env = SplitBenders.CFLPSplitSubEnv(data,algo_params)
check_df = DataFrame(Instance = String[], Check = Bool[])
# for i in 1:66
# for i in 67:71
for i in 30:30
    instance = "p$i"
    @info "Instance: $instance"
    data = SplitBenders.read_data(instance)

    master_env = SplitBenders.MasterProblem(data, solver=solver)
    relax_integrality(master_env.model)
    sub_env = SplitBenders.CFLPSplitSubEnv(data,algo_params, solver=solver) 
    df = SplitBenders.run_Benders(data,master_env,sub_env)
    obj_split = df[end,:LB]

    mip_env = SplitBenders.CFLPMipEnv(data)
    optimize!(mip_env.model)
    obj_mip = objective_value(mip_env.model)

    check = obj_split ≈ obj_mip
    push!(check_df, (instance, check))
    if !check @error "Instance: $instance: Mismatch in objective values: $obj_split vs $obj_mip" end
end

# CSV.write("test/check.csv", check_df)




