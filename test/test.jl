include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV, CPLEX, Gurobi

# instance = "f50-c50-r5.0-p1"
# data = DisjunctiveBenders.read_random_data(instance)

instance = "ga250a-5"
data = SplitBenders.read_orlib_file(instance)
@info data
