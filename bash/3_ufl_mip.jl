include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV, Logging, DataFrames


#-----------------------------------------------------------------------
solver = :Gurobi
# solver = :CPLEX

settings = SplitBenders.parse_commandline()
instance = settings["instance"]
data = SplitBenders.read_Simple_data(instance; filepath = "src/BendersDatasets/KoerkelGhosh-asym/")


df = DataFrame(LB = Float64[], UB = Float64[], gap = Float64[])
mip_env = SplitBenders.UFLPMipEnv(data)
set_time_limit_sec(mip_env.model, 7200)
# optimize!(mip_env.model)

io = open("results4/MIP/result_$(instance).txt", "w+")
logger = SimpleLogger(io)
with_logger(logger) do
    optimize!(mip_env.model)
    LB = JuMP.objective_bound(mip_env.model)
    UB = JuMP.objective_value(mip_env.model)
    gap = JuMP.relative_gap(mip_env.model)
    new_row = (LB,UB,gap)
    push!(df,new_row)
    CSV.write("results4/MIP/result_$(instance).csv", df)
end
flush(io)
close(io)
# df = SplitBenders.run_Benders(data,master_env,sub_env)

# result post processing
# CSV.write("temp/result_$(instance)_$(cut_strategy)_$(SplitCGLPNormType)_$(SplitSetSelectionPolicy)_$(StrengthenCutStrategy)_$(SplitBendersStrategy)_2.csv", df)





#--------------------------------------------------------------
# cut_strategy = settings["cut_strategy"]
# SplitCGLPNormType = settings["SplitCGLPNormType"]
# SplitSetSelectionPolicy = settings["SplitSetSelectionPolicy"]
# StrengthenCutStrategy = settings["StrengthenCutStrategy"]
# SplitBendersStrategy = settings["SplitBendersStrategy"]

# algo_params = SplitBenders.AlgorithmParams(
#     cut_strategy = cut_strategy,
#     SplitCGLPNormType = SplitCGLPNormType,
#     SplitSetSelectionPolicy = SplitSetSelectionPolicy,
#     StrengthenCutStrategy = StrengthenCutStrategy,
#     SplitBendersStrategy = SplitBendersStrategy
# )

# master_env = SplitBenders.MasterProblem(data)
# sub_env = SplitBenders.CFLPSplitSubEnv(data,algo_params)

# df = SplitBenders.run_Benders(data,master_env,sub_env)

# # result post processing
# cut_strategy = SplitBenders.reConvert_cutstrategy(cut_strategy)
# SplitCGLPNormType = SplitBenders.reConvert_normtype(SplitCGLPNormType)
# SplitSetSelectionPolicy = SplitBenders.reConvert_splitsetselectionpolicy(SplitSetSelectionPolicy)
# StrengthenCutStrategy = SplitBenders.reConvert_splitstengtheningpolicy(StrengthenCutStrategy)
# SplitBendersStrategy = SplitBenders.reConvert_splitbenderspolicy(SplitBendersStrategy)
# CSV.write("results/Gurobi/result_$(instance)_$(cut_strategy)_$(SplitCGLPNormType)_$(SplitSetSelectionPolicy)_$(StrengthenCutStrategy)_$(SplitBendersStrategy).csv", df)
