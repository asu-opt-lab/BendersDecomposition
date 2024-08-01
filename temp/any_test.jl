include("../src/SplitBenders.jl")
import .SplitBenders
using JuMP, CSV, Logging, DataFrames, CPLEX, Gurobi

instance = "MO1"
@info "Instance: $instance"
data = SplitBenders.read_Orlib_data(instance; filepath = "src/BendersDatasets/M/O")