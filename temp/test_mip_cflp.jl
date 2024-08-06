include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV, Logging, DataFrames, CPLEX, Gurobi

solver = :Gurobi
# solver = :CPLEX

instance = "f700-c700-r5.0-p1"
data = SplitBenders.read_GK_data(instance)


mip_env = SplitBenders.CFLPMipEnv(data)
optimize!(mip_env.model)