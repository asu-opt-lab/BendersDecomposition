include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV, Logging, DataFrames, CPLEX, Gurobi

solver = :Gurobi

instance = "ga250a-3"
@info "Instance: $instance"
data = SplitBenders.read_Simple_data(instance; filepath = "src/BendersDatasets")

mip_env = SplitBenders.UFLPMipEnv(data)
optimize!(mip_env.model)