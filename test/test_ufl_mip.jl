include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV, CPLEX, Gurobi

# instance = "f50-c50-r5.0-p1"
# data = DisjunctiveBenders.read_random_data(instance)

instance = "B1.3"
data = SplitBenders.read_data_UFL(instance)
@info instance


mip_env = SplitBenders.UFLPMipEnv(data)
optimize!(mip_env.model)
@show objective_value(mip_env.model)