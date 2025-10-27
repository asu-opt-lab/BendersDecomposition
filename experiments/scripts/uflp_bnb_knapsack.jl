using JuMP, DataFrames, Logging, CSV
using BendersDecomposition
using Random
using Printf  
using Statistics  
import BendersDecomposition: generate_cuts
include("$(dirname(@__DIR__))/example/uflp/data_reader.jl")
include("$(dirname(@__DIR__))/example/uflp/oracle.jl")
include("$(dirname(@__DIR__))/example/uflp/model.jl")

# load settings
args = parse_commandline()

Random.seed!(args["seed"])
instance = args["instance"]
output_dir = args["output_dir"]

# -----------------------------------------------------------------------------
# load problem data
# -----------------------------------------------------------------------------
problem = read_Simple_data(instance)
dim_x = problem.n_facilities
dim_t = problem.n_customers
c_x = problem.fixed_costs
# c_x = zeros(dim_x)
c_t = ones(dim_t)
data = Data(dim_x, dim_t, problem, c_x, c_t)

# -----------------------------------------------------------------------------
# load parameters
# -----------------------------------------------------------------------------
# Algorithm parameters
benders_param = BendersBnBParam(
    time_limit = 14400.0,
    gap_tolerance = 1e-6,
    #disjunctive_root_process = true,
    verbose = true
)

# Solver parameters
master_solver_param = Dict(
    "solver" => "CPLEX", 
    "CPX_PARAM_EPINT" => 1e-9, 
    "CPX_PARAM_EPRHS" => 1e-9,
    "CPX_PARAM_EPGAP" => 1e-6,
    "CPX_PARAM_BRDIR" => 1
)

# -----------------------------------------------------------------------------
# master model
# -----------------------------------------------------------------------------
master = Master(data; solver_param = master_solver_param)
update_model!(master, data)

# -----------------------------------------------------------------------------
# typical oracles
# -----------------------------------------------------------------------------
# Create two oracles for kappa & nu
typical_oracle = UFLKnapsackOracle(data)

for k=1:2
    set_parameter!(typical_oracle, "add_only_violated_cuts", true)
end

# -----------------------------------------------------------------------------
# root node preprocessing
# -----------------------------------------------------------------------------
root_seq_type = BendersSeq
root_param = BendersSeqParam(
    time_limit = 100.0,
    gap_tolerance = 1e-9,
    verbose = true
)

# Create root node preprocessing with oracle
root_preprocessing = RootNodePreprocessing(typical_oracle, root_seq_type, root_param)

# -----------------------------------------------------------------------------
# lazy callback
# -----------------------------------------------------------------------------
lazy_callback = LazyCallback(typical_oracle)

# -----------------------------------------------------------------------------
# user callback
# -----------------------------------------------------------------------------
user_callback = NoUserCallback()

# -----------------------------------------------------------------------------
# BendersBnB
# -----------------------------------------------------------------------------
env = BendersBnB(
    data, 
    master, 
    root_preprocessing, 
    lazy_callback, 
    user_callback; 
    param = benders_param
)

# -----------------------------------------------------------------------------
# solve
# -----------------------------------------------------------------------------
solution_log = solve!(env)





